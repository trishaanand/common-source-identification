/* Use assertions. */
use Assert;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

/* For ceiling function */
use Math;

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
config const numThreads : int = 10; 

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
  const numDomain : domain(1) = {1..n};
  const threadDomain : domain(1) = {0..#numThreads};

  var corrMatrix : [corrDomain] real;
  var nrCorrelations = (n * (n - 1)) / 2;
  var crossDomain : domain(1) = {0..#nrCorrelations};
  // Map to maintain relation b/w index in corrMatrix and images required to calculate pce for that index.
  var crossTuples : [crossDomain] 2*int;

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

  // Create a timer for each thread. Since all the threads will execute parallely, 
  // we'll get the max of these timers to get the overall time.
  var threadTimer : [threadDomain] Timer;

  class ThreadData {
    var prnu, prnuRot, resultComplex : [imageDomain] complex;
    var fwPlan : fftw_plan;
    var fwPlanRot : fftw_plan;
    var bwPlan : fftw_plan;

    proc init(h : int, w : int, id : int, ref data : prnu_data) {
      prnuInit(h, w, data);
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

  var threadArray : [threadDomain] ThreadData;
  var data : [threadDomain] prnu_data;

  /* This is the number of correlations that must be performed by each thread
   * 1. Divide the total # of calculations by the number of threads. 
   * 2. If remainder is non-zero, distribute the remaining calculations 1 per thread
   * This ensures efficient work distribution of computations amongst all the threads.
   */
  var defaultNum = nrCorrelations/numThreads : int;
  var rem = nrCorrelations % numThreads;
  var high : int;
  var threadTuples : [threadDomain] 2*int;
  high = -1; 

  /* IMP: Must be sequential else it throws a segmentation fault.
   * This is because the FFT & PRNU  initialization must be performed sequentially
   */
  for thread in threadDomain {
    // Calculate the indices of the the corrMatrix that each thread will compute over
    var low = high + 1; 
    var num = defaultNum;
    if thread < rem {
      num = defaultNum + 1;
    }
    high = low + num -1;
    threadTuples(thread) = (low, high);
    // Init all the data for each thread.
    threadArray[thread] = new unmanaged ThreadData(h, w, thread, data[thread]);
  }

  coforall thread in threadDomain {
    // Start the thread timer
    threadTimer[thread].start();

    var prnuTemp : [imageDomain] complex;
    var images : [2][imageDomain] RGB;

    var localSubDom : domain(1) = {threadTuples(thread)[1]..threadTuples(thread)[2]};
    threadTimer[thread].stop();

    // Sequentially iterate over the localSubdom. MUST BE SEQUENTIAL because of FFT crap
    for idx in localSubDom {
      var (i,j) = crossTuples(idx);

      // Read both the images from disk 
      readJPG(images[0], imageFileNames[i]);
      readJPG(images[1], imageFileNames[j]);

      // Start the thread timer
      threadTimer[thread].start();

      calculatePrnuComplex(h, w, images[0], threadArray[thread].prnu, data[thread]);
      calculatePrnuComplex(h, w, images[1], prnuTemp, data[thread]);
      rotated180Prnu(h, w, prnuTemp, threadArray[thread].prnuRot);

      // Calculate the FFT on the prnu & rot arrays
      execute(threadArray[thread].fwPlan);
      execute(threadArray[thread].fwPlanRot);

      // Calculate the dot product
      threadArray[thread].resultComplex = threadArray[thread].prnu * threadArray[thread].prnuRot;
      
      // Inverse FFT
      execute(threadArray[thread].bwPlan);

      // Compute the pce value
      corrMatrix(i,j) = computePCE(h, w, threadArray[thread].resultComplex);

      threadTimer[thread].stop();
    }
  }

  cleanup();
  
  var overallTimer = max reduce threadTimer.elapsed();

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer, "s");

  writeln("Throughput: ", nrCorrelations / overallTimer, " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}
