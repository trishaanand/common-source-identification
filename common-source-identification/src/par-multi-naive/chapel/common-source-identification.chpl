/* Use assertions. */
use Assert;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

use BlockDist;

/* FFT Module */
use FFTW_MT;

/* Read in JPG images. */
use read_jpg;

/* Compute PRNU noise patterns. */
use prnu;

/* user defined functions */
use pce;
use utils;

/* Configuration parameters */
config const imagedir : string = "images";
config const writeOutput : bool = false;
config const numThreads : int = 8; 

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnuComplex(h : int, w : int, image : [] RGB, prnuComplex : [] complex, ref data : prnu_data) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};

  var prnu : [imageDomain] real;
  
  prnuExecute(prnu, image, data);

  forall (i,j) in imageDomain {
    prnuComplex(i,j) = prnu(i,j) + 0i;
  }
}

/* Rotate the prnu matrix by 180 degrees */
proc rotated180Prnu(h : int, w : int, prnu : [] complex, prnuRot : [] complex) {
  const imageDomain: domain(2) = {0..#h, 0..#w};

  forall (i,j) in imageDomain do 
    prnuRot(i,j) = prnu(h-i-1, w-j-1);
}

proc main() {
  run();
}

proc run() {
  /* Obtain the images. */
  var imageFileNames = getImageFileNames(imagedir);

  /* n represents the number of images that have to be correlated. */
  var n = imageFileNames.size;

  /* h, w will represent the height and width of an image or PRNU noise pattern 
   * throughout the code.
   */
  var h, w : int;
  (h, w) = getDimensionsJPG(imageFileNames.front());

  /* Create a domain for the correlation matrix. */
  const corrDomain : domain(2) = {1..n, 1..n};
  const imageDomain : domain(2) = {0..#h, 0..#w};
  const threadDomain : domain(1) = {0..#numThreads};

  var corrMatrix : [corrDomain] real;
  var nrCorrelations = (n * (n - 1)) / 2;
  var crossDomain = {0..#nrCorrelations} dmapped Block({0..#nrCorrelations});
  var crossTuples : [crossDomain] 2*int;
  const localeDomain = {0..#numLocales} dmapped Block ({0..#numLocales});
  var overallTimerLoc : [localeDomain] real;
  
  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  writeln("  ", numThreads, " numThreads");

  /* ************************* Start here ********************* */
  forall (i,j) in corrDomain {
    if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
    }
  }

  class ThreadData {
    var prnu, prnuRot, resultComplex : [imageDomain] complex;
    var fwPlan : fftw_plan;
    var fwPlanRot : fftw_plan;
    var bwPlan : fftw_plan;

    proc init(h_loc : int, w_loc : int, id : int, ref data : prnu_data) {
      prnuInit(h_loc, w_loc, data);
      
      fwPlan = plan_dft(prnu, prnu, FFTW_FORWARD, FFTW_ESTIMATE);
      fwPlanRot = plan_dft(prnuRot, prnuRot, FFTW_FORWARD, FFTW_ESTIMATE);
      bwPlan = plan_dft(resultComplex, resultComplex, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    proc deinit() {
      destroy_plan(fwPlan);
      destroy_plan(fwPlanRot);
      destroy_plan(bwPlan);
    }
  }

  coforall loc in Locales do on loc {
    // Create a timer for each thread. We'll sum these timers to get the overall time
    var localeTimer : Timer;
    var crossSubDomain = crossDomain.localSubdomain();
    var threadArray : [threadDomain] ThreadData;
    var data : [threadDomain] prnu_data;

    /* Calculate the indices of the the corrMatrix that each thread will compute over
     * We attempt to distribute the computations evenly over all the threads. 
     * Any remaining computations are then distributed sequentially to threads to achieve optimal
     * computation distribution.
     */
    var defaultNum = crossSubDomain.size/numThreads : int;
    var rem = crossSubDomain.size % numThreads;
    var high : int;
    var threadTuples : [threadDomain] 2*int;
    high = crossSubDomain.low - 1; 

    // Calculate the indices of the the corrMatrix that each thread will compute over
    // IMP: Must be sequential else it throws a segmentation fault
    for thread in threadDomain {
      var low = high + 1; 
      var num = defaultNum;
      if thread < rem {
        num = defaultNum + 1;
      }
      high = low + num -1;
      threadTuples(thread) = (low, high);
      // Create unmanaged objects for each thread. Data is initialized in the init fxn of the class.
      threadArray[thread] = new unmanaged ThreadData(h, w, thread, data[thread]);
    }

    // Start the locale timer
    localeTimer.start();
    coforall thread in threadDomain {

      var prnuTemp : [imageDomain] complex;
      var localSubDom : domain(1) = {threadTuples(thread)[1]..threadTuples(thread)[2]};
      
      // Sequentially iterate over the localSubdom. MUST BE SEQUENTIAL because of FFT module limitations
      for idx in localSubDom {
        var (i,j) = crossTuples(idx);

        // Read both the images from disk 
        var image, imageRot : [imageDomain] RGB;
        readJPG(image, imageFileNames[i].localize());
        readJPG(imageRot, imageFileNames[j].localize());

        // Calculate the prnu for both images and rotate the prnu for 2nd image
        calculatePrnuComplex(h, w, image, threadArray[thread].prnu, data[thread]);
        calculatePrnuComplex(h, w, imageRot, prnuTemp, data[thread]);
        rotated180Prnu(h, w, prnuTemp, threadArray[thread].prnuRot);

        // Calculate the FFT on the prnu & rot arrays
        execute(threadArray[thread].fwPlan);
        execute(threadArray[thread].fwPlanRot);

        // Calculate the dot product
        threadArray[thread].resultComplex = threadArray[thread].prnu * threadArray[thread].prnuRot;
        
        // Inverse FFT
        execute(threadArray[thread].bwPlan);

        // Compute the PCE value
        corrMatrix(i,j) = computePCE(h, w, threadArray[thread].resultComplex);
      }
    }
    localeTimer.stop();
    overallTimerLoc[loc.id] = localeTimer.elapsed();

    for thread in threadDomain {
      delete threadArray[thread];
    }
    cleanup();
  }
  
  var overallTimer = max reduce overallTimerLoc;
  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer, "s");

  writeln("Throughput: ", nrCorrelations / overallTimer, " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}
