/* Use assertions. */
use Assert;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

use BlockDist;

use VisualDebug;

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
config const maxCache : int = 5;

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnuComplex(h : int, w : int, image : [] RGB, prnuComplex : [] complex, ref data : prnu_data) {
  
  /* Create a domain for an image and allocate the image itself */
  var imageDomain: domain(2) = {0..#h, 0..#w};

  var prnu : [imageDomain] real;
  
  prnuExecute(prnu, image, data);

  prnuComplex = prnu;
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
  flushWriteln("  ", maxCache, " maxCache");

  /* ************************* Start here ********************* */
  forall (i,j) in corrDomain {
    if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
    }
  }

  class ThreadData {
    var prnu, prnuRot, resultComplex : [imageDomain] complex;
    var i, j : int;
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
    // tagVdebug("Coforall Locales");
    var crossSubDomain = crossDomain.localSubdomain();
    var localNumThreads = numThreads;
    if (crossSubDomain.size < localNumThreads) {
      localNumThreads = crossSubDomain.size;
    }
    var localeTimer : Timer;
    if (localNumThreads != 0) {
      // This is the number of correlations that must be performed by each thread
      var threadDomain : domain(1) = {0..#localNumThreads};
      // Create a timer for each thread. We'll sum these timers to get the overall time
      var threadArray : [threadDomain] ThreadData;
      var data : [threadDomain] prnu_data;
      var h_loc, w_loc : int;
      h_loc = h;
      w_loc = w;
      var imageDomainLoc : domain(2) = {0..#h_loc, 0..#w_loc};
      var defaultNum = crossSubDomain.size/localNumThreads : int;
      var rem = crossSubDomain.size % localNumThreads;
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
        prnuInit(h_loc, w_loc, data[thread]);
        threadArray[thread] = new unmanaged ThreadData();
        threadArray[thread].fwPlan = plan_dft(threadArray[thread].prnu, threadArray[thread].prnu, FFTW_FORWARD, FFTW_ESTIMATE);
        threadArray[thread].fwPlanRot = plan_dft(threadArray[thread].prnuRot, threadArray[thread].prnuRot, FFTW_FORWARD, FFTW_ESTIMATE);
        threadArray[thread].bwPlan = plan_dft(threadArray[thread].resultComplex, threadArray[thread].resultComplex, FFTW_BACKWARD, FFTW_ESTIMATE);
      }

      /* Create a cache for prnu & prnuRot 
       * This is implemented via a arrays over the maxCache of images that are required
       * in the correlation. After each calculation, we'll save the prnu & prnuRot in the cache. 
       * If the max cap (arbitrary number) of the cache is reached, then don't save it. 
       * This is to save on memory.
       */
      var cachePrnuIdx, cachePrnuRotIdx  : [{0..#maxCache}] int;
      var cachePrnu : [{0..#maxCache}][imageDomain] complex;
      var cachePrnuRot : [{0..#maxCache}][imageDomain] complex;
      var cachePrnuSize, cachePrnuRotSize : atomic int;
      var counterPrnu, counterPrnuRot : atomic int;
      localeTimer.start();
      coforall thread in threadDomain {
        // Start the thread timer
        // threadTimer[thread].start();
        var prnuTemp : [imageDomainLoc] complex;

        var localSubDom : domain(1) = {threadTuples(thread)[1]..threadTuples(thread)[2]};
        // threadTimer[thread].stop();
        // Sequentially iterate over the localSubdom. MUST BE SEQUENTIAL because of FFT crap
        for idx in localSubDom {
          var (i,j) = crossTuples(idx);

          var flagI, flagJ : bool;
          var flagIdx, flagJdx : int;
          var image, imageRot : [imageDomain] RGB;

          (flagI, flagIdx) = cachePrnuIdx.find(i);
          (flagJ, flagJdx) = cachePrnuRotIdx.find(j);

          //Check for cache hit
          if(flagI) {
            // flushWriteln("On locale: ", here.id, " Read prnu from cache for i: ", i, " flagIdx: ", flagIdx);
            threadArray[thread].prnu = cachePrnu[flagIdx];
            counterPrnu.add(1);
          } else {
            readJPG(image, imageFileNames[i].localize());
            calculatePrnuComplex(h_loc, w_loc, image, threadArray[thread].prnu, data[thread]);
            //  Calculate the FFT on the prnu arrays
            execute(threadArray[thread].fwPlan);
          }
          
          //Check for cache hit
          if(flagJ) {
            // flushWriteln("On locale: ", here.id , " Read prnuRot from cache for j: ", j , " flagJdx: ", flagJdx);
            threadArray[thread].prnuRot = cachePrnuRot[flagJdx];
            counterPrnuRot.add(1);
          } else {
            readJPG(imageRot, imageFileNames[j].localize());
            calculatePrnuComplex(h_loc, w_loc, imageRot, prnuTemp, data[thread]);
            rotated180Prnu(h_loc, w_loc, prnuTemp, threadArray[thread].prnuRot);
            //  Calculate the FFT on rot arrays
            execute(threadArray[thread].fwPlanRot);
          }

          // Calculate the dot product
          threadArray[thread].resultComplex = threadArray[thread].prnu * threadArray[thread].prnuRot;
          
          // Inverse FFT
          execute(threadArray[thread].bwPlan);

          corrMatrix(i,j) = computePCE(h_loc, w_loc, threadArray[thread].resultComplex);
          threadArray[thread].i = i;
          threadArray[thread].j = j;

          // Write the values to cache so that we don't need to calculate it again
          if (cachePrnuSize.read() < (maxCache - 1)) {
            var (found, val) = cachePrnuIdx.find(i);
            // flushWriteln("On locale: ", here.id, " In the prnu cache, got found: ", found, " and val: ", val, " for i: ", i);
            if (!found) {
              cachePrnuIdx[cachePrnuSize.read()] = i;
              cachePrnu[cachePrnuSize.read()] = threadArray[thread].prnu;
              cachePrnuSize.add(1);
              // flushWriteln("On locale: ", here.id, " Writing i: ", i, " cachePrnuSize: ", cachePrnuSize );
            }
          }
          if (cachePrnuRotSize.read() < (maxCache - 1)) {
            var (found, val) = cachePrnuRotIdx.find(j);
            // flushWriteln("On locale: ", here.id, " In the prnuRot cache, got found: ", found, " and val: ", val, " for j: ", j);
            if (!found) {
              cachePrnuRotIdx[cachePrnuRotSize.read()] = j;
              cachePrnuRot[cachePrnuRotSize.read()] = threadArray[thread].prnuRot;
              cachePrnuRotSize.add(1);
              // flushWriteln("On locale: ", here.id, " Writing j: ", j, " cachePrnuRotSize: ", cachePrnuRotSize);
            }
          }
          // threadTimer[thread].stop();
        }
        // threadTimer[thread].stop();
        // Cleanup everything
        // delete threadArray[thread];
      }
      localeTimer.stop();
      overallTimerLoc[loc.id] = localeTimer.elapsed();
      flushWriteln("On locale: ", here.id, " cachePrnuSize: ", cachePrnuSize, ", cache hits : ", counterPrnu );
      flushWriteln("On locale: ", here.id, " cachePrnuRotSize: ", cachePrnuRotSize, ", cache hits : ", counterPrnuRot );
      // for thread in threadDomain {
      //   delete threadArray[thread];
      // }
      // cleanup();
    }
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
