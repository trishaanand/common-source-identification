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

/* user defined modules */
use fft;
use utils;

/* Configuration parameters */
config const imagedir : string = "images";
config const writeOutput : bool = false;
config const numThreads : int = 8; 
config const maxCache : int = 10;

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnu(imageDomain : domain, image : [] RGB, prnuComplex : [] complex, ref data : prnu_data) {
  
  var prnu : [imageDomain] real;
  prnuExecute(prnu, image, data);
  prnuComplex = prnu;
}

/* Rotate the matrix 180 degrees */
proc rotatePrnuBy180(h : int, w : int, prnu : [] complex, prnuRot : [] complex) {
  const imageDomain: domain(2) = {0..#h, 0..#w};

  forall (i,j) in imageDomain do 
    prnuRot(h-i-1, w-j-1) = prnu(i,j);
}

proc writeCache(cacheSize : atomic int, i : int, ref cacheIdx : [] int, ref cache, ref val : [] complex) {
    if (cacheSize.read() < (maxCache - 1)) {
    var (found, idxVal) = cacheIdx.find(i);
    if (!found) {
      cacheIdx[cacheSize.read()] = i;
      cache[cacheSize.read()] = val;
      cacheSize.add(1);
    }
  }
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

  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  writeln("  ", numThreads, " numThreads");
  writeln("  ", maxCache, " maxCache");

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
    var localeTimer : Timer;

    /* localNumThreads determines the number of threads on each locale. 
     * If localNumThreads == 0, that means there are no threads required to be spawned on that locale. 
     * Required if we have more threads than number of correlations to compute.
     */    
    var crossSubDomain = crossDomain.localSubdomain();
    var localNumThreads = numThreads;
    if (crossSubDomain.size < localNumThreads) {
      localNumThreads = crossSubDomain.size;
    }
    
    if (localNumThreads != 0) {

      // This is the number of correlations that must be performed by each thread
      var threadDomain : domain(1) = {0..#localNumThreads};

      var threadArray : [threadDomain] ThreadData;
      var data : [threadDomain] prnu_data;
      var h_loc, w_loc : int;
      h_loc = h;
      w_loc = w;
      var imageDomainLoc : domain(2) = {0..#h_loc, 0..#w_loc};

      /* Calculate the indices of the the corrMatrix that each thread will compute over
       * IMP: Must be sequential else it throws a segmentation fault
       */
      var defaultNum = crossSubDomain.size/localNumThreads : int;
      var rem = crossSubDomain.size % localNumThreads;
      var high : int;
      var threadTuples : [threadDomain] 2*int;
      high = crossSubDomain.low - 1; 

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
      
      // Start the timer for this locale
      localeTimer.start();
      
      coforall thread in threadDomain {
        var localSubDom : domain(1) = {threadTuples(thread)[1]..threadTuples(thread)[2]};
        
        /* Sequentially iterate over the localSubdom. MUST BE SEQUENTIAL because of how FFT plans works.  
         * FFTs should be calculated on the same thread as the plan
         */
        for idx in localSubDom {
          var (i,j) = crossTuples(idx);

          var flagI, flagJ : bool;
          var flagIdx, flagJdx : int;
          var image, imageRot : [imageDomain] RGB;

          (flagI, flagIdx) = cachePrnuIdx.find(i);
          (flagJ, flagJdx) = cachePrnuRotIdx.find(j);

          /* If the prnu value is not found in cache, then calculate it */
          if(!flagI) {
            readJPG(image, imageFileNames[i].localize());
            calculatePrnu(imageDomainLoc, image, threadArray[thread].prnu, data[thread]);
            //  Calculate the FFT on the prnu arrays
            execute(threadArray[thread].fwPlan);
          }

          /* If the prnuRot value is not found in cache, then calculate it */
          if(!flagJ) {
            var prnuTemp : [imageDomainLoc] complex;
            readJPG(imageRot, imageFileNames[j].localize());
            calculatePrnu(imageDomainLoc, imageRot, prnuTemp, data[thread]);
            rotatePrnuBy180(h_loc, w_loc, prnuTemp, threadArray[thread].prnuRot);
            //  Calculate the FFT on rot arrays
            execute(threadArray[thread].fwPlanRot);
          }

          /* We use pointers to the prnu & prnuRot array values so that we don't have to copy the actual 
           * values over from cache to the threadArray class object.
           */
          ref prnuRef = if (flagI) then cachePrnu[flagIdx] else threadArray[thread].prnu;
          ref prnuRotRef = if (flagJ) then cachePrnuRot[flagJdx] else threadArray[thread].prnuRot;

          // Calculate the dot product
          forall (x,y) in imageDomainLoc {
            threadArray[thread].resultComplex(x,y) = prnuRef(x,y) * prnuRotRef(x,y);
          }
          
          // Inverse FFT
          execute(threadArray[thread].bwPlan);

          computePCE(h_loc, w_loc, threadArray[thread].resultComplex, corrMatrix(i,j));

          // Write the values to cache so that we don't need to calculate it again
          writeCache(cachePrnuSize, i, cachePrnuIdx, cachePrnu, prnuRef);
          writeCache(cachePrnuRotSize, j, cachePrnuRotIdx, cachePrnuRot, prnuRotRef);
        }
      }
      localeTimer.stop();
      overallTimerLoc[loc.id] = localeTimer.elapsed();
    }
  }

  var overallTimer = max reduce overallTimerLoc;
  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer, "s");

  writeln("Throughput: ", nrCorrelations / overallTimer, " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }
}
