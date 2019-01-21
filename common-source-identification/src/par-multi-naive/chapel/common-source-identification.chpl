/* Use assertions. */
use Assert;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

use BlockDist;

/* Read in JPG images. */
use read_jpg;

/* Compute PRNU noise patterns. */
use prnu;

/* user defined functions */
use fft;
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

proc rotated180Prnu(h : int, w : int, prnu : [] complex, prnuRot : [] complex) {
  const imageDomain: domain(2) = {0..#h, 0..#w};

  /* Rotate the matrix 180 degrees */
  forall (i,j) in imageDomain do 
    prnuRot(i,j) = prnu(h-i-1, w-j-1);
}

proc readImage(i: int, imageDomain:domain) {
  var image : [imageDomain]RGB;
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
  const numDomain : domain(1) = {1..n};
  const threadDomain : domain(1) = {0..#numThreads};

  var corrMatrix : [corrDomain] real;
  var nrCorrelations = (n * (n - 1)) / 2;
  var crossDomain = {0..#nrCorrelations} dmapped Block({0..#nrCorrelations});
  var crossTuples : [crossDomain] 2*int;
  const localeDomain = {0..#numLocales} dmapped Block ({0..#numLocales});
  var overallTimerLoc : [localeDomain] real;
  
  flushWriteln("Running Common Source Identification...");
  flushWriteln("  ", n, " images");
  flushWriteln("  ", numLocales, " locale(s)");
  flushWriteln("  ", numThreads, " numThreads");

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

    proc deinit() {
      // prnuDestroy(data);
      destroy_plan(fwPlan);
      destroy_plan(fwPlanRot);
      destroy_plan(bwPlan);
    }
  }

  coforall loc in Locales do on loc {
    // Create a timer for each thread. We'll sum these timers to get the overall time
    var threadTimer : [threadDomain] Timer;
    var crossSubDomain = crossDomain.localSubdomain();
    var threadArray : [threadDomain] ThreadData;
    var data : [threadDomain] prnu_data;

    // This is the number of correlations that must be performed by each thread
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
      // Init all the data
      prnuInit(h, w, data[thread]);
      threadArray[thread] = new unmanaged ThreadData();
      threadArray[thread].fwPlan = plan_dft(threadArray[thread].prnu, threadArray[thread].prnu, FFTW_FORWARD, FFTW_ESTIMATE);
      threadArray[thread].fwPlanRot = plan_dft(threadArray[thread].prnuRot, threadArray[thread].prnuRot, FFTW_FORWARD, FFTW_ESTIMATE);
      threadArray[thread].bwPlan = plan_dft(threadArray[thread].resultComplex, threadArray[thread].resultComplex, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    coforall thread in threadDomain {
      // Start the thread timer
      threadTimer[thread].start();

      var prnuTemp : [imageDomain] complex;

      var localSubDom : domain(1) = {threadTuples(thread)[1]..threadTuples(thread)[2]};
      threadTimer[thread].stop();
      // Sequentially iterate over the localSubdom. MUST BE SEQUENTIAL because of FFT crap
      for idx in localSubDom {
        var (i,j) = crossTuples(idx);

        // Read both the images from disk 
        var image, imageRot : [imageDomain] RGB;
        readJPG(image, imageFileNames[i].localize());
        readJPG(imageRot, imageFileNames[j].localize());

        // Start the thread timer
        threadTimer[thread].start();

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

        corrMatrix(i,j) = computePCE(h, w, threadArray[thread].resultComplex);
        threadTimer[thread].stop();
      }
      // Cleanup everything
      // delete threadArray[thread];
    }
    overallTimerLoc[loc.id] = max reduce threadTimer.elapsed();

    for thread in threadDomain {
      delete threadArray[thread];
    }
    cleanup();
  }
  
  var overallTimer = max reduce overallTimerLoc;
  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer, "s");

  flushWriteln("Throughput: ", nrCorrelations / overallTimer, " corrs/s");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}
