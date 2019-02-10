/* Use assertions. */
use Assert;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

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
config const num : int = 10; 

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnuComplex(h : int, w : int, image : [] RGB, ref data : prnu_data, prnuComplex : [] complex) {
  
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

proc main() {
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
  var corrMatrix : [corrDomain] real;
  var nrCorrelations = (n * (n - 1)) / 2;
  var crossDomain : domain(1) = {0..#nrCorrelations};

  var crossTuples : [crossDomain] 2*int;
  var overallTimer, prnuTimer, fftTimer, corrTimer, crossTimer : Timer;

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const numDomain : domain(1) = {1..n};

  flushWriteln("Running Common Source Identification...");
  flushWriteln("  ", n, " images");
  flushWriteln("  ", numLocales, " locale(s)");

  /* ************************* Start here ********************* */
  //crossTuple stores the mapping between the cross domain and images required.
  //i is for prnu, j would be for rotation.
  forall (i,j) in corrDomain {
    if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
    }
  }

  var idx : int = 0;
  
  var data : [numDomain] prnu_data;  
  for i in numDomain {
    prnuInit(h, w, data(i));
  }
  printGlobalMemory("******* After prnuInit ******");
  //Sleep is used to manually read /proc/meminfo on the locale.
  // sleep(60);

  while(idx < nrCorrelations) {
    printGlobalMemory("******* Beginning of while loop ******");
    // sleep(60);
    // printGlobalMemory("******* Continuing after sleep in while loop ******  ");

    overallTimer.start();

    //Calculating the range of correlations for this batch iteration.
    var localRange : int;
    if (idx + num) >= nrCorrelations {
      localRange = nrCorrelations-1;
    } else {
      localRange = idx + num - 1;
    }

    const localSubDom : domain(1) = {idx..localRange};
    //Sparse domains store the indices required for images, prnus and rotated prnus for calculating all
    //the correlations for this batch.
    var imgSparseDom, prnuSparseDom, prnuRotSparseDom : sparse subdomain(numDomain);
    
    // Sequential iteration because we can't add indices to sparse in forall
    for it in localSubDom {
      var (i,j) = crossTuples(it);
      //Append to the sparse domain index
      imgSparseDom += i;
      imgSparseDom += j;
      prnuSparseDom += i; 
      prnuRotSparseDom += j;
    }

    /* Perform all initializations here. Hence stopping the timer */
    overallTimer.stop();
    var images : [imgSparseDom][imageDomain] RGB;
    var prnuArray : [prnuSparseDom][imageDomain] complex;
    var prnuRotArray : [prnuRotSparseDom][imageDomain] complex;
    var fftPlan : [prnuSparseDom] fftw_plan;
    var fftRotPlan : [prnuRotSparseDom] fftw_plan;

    for i in imgSparseDom {
      readJPG(images[i], imageFileNames[i]);
    }

    /* Start the timer here */
    overallTimer.start();

    prnuTimer.start();
    //for all images required in this batch, calculate prnu.
    forall i in imgSparseDom {
      var prnu : [imageDomain] complex;
      calculatePrnuComplex(h, w, images[i], data(i), prnu);
      if prnuSparseDom.member(i) {
        prnuArray(i) = prnu;
      }
      if prnuRotSparseDom.member(i) {
        rotated180Prnu(h, w, prnu, prnuRotArray(i));
      }
    }
    prnuTimer.stop();

    fftTimer.start();
    for i in prnuSparseDom {
      fftPlan(i) = planFFT(prnuArray(i), FFTW_FORWARD);  
    }
    for j in prnuRotSparseDom {
      fftRotPlan(j) = planFFT(prnuRotArray(j), FFTW_FORWARD);
    }

    sync {
      begin {
        forall i in prnuSparseDom {
          execute(fftPlan(i));
        }
      }
      begin {
        forall i in prnuRotSparseDom {
          execute(fftRotPlan(i));
        }
      }
    }
    fftTimer.stop();
    //Destroy the fft plans
    for i in prnuSparseDom {
      destroy_plan(fftPlan(i));
    }
    for i in prnuRotSparseDom {
      destroy_plan(fftRotPlan(i));
    }

    crossTimer.start();
    var resultComplex : [localSubDom][imageDomain] complex;

    //Calculate dot product of prnu and rotated prnu array and save it in resultComplex
    forall it in localSubDom {
      var (i,j) = crossTuples(it);
      forall (x, y) in imageDomain {
        resultComplex(it)(x,y) = prnuArray(i)(x,y) * prnuRotArray(j)(x,y);
      }
    } 
    crossTimer.stop();

    corrTimer.start();
    var fftPlanBack : [localSubDom] fftw_plan;

    // Plan the inverse fft
    for it in localSubDom {
      fftPlanBack(it) = planFFT(resultComplex(it), FFTW_BACKWARD) ; 
    }

    forall it in localSubDom {
      execute(fftPlanBack(it));
    }

   for it in localSubDom {
      destroy_plan(fftPlanBack(it));
    }

    /* scale and compute pce now */
    forall it in localSubDom {
      var (i,j) = crossTuples(it);
      
      var result : [imageDomain] real = (resultComplex(it).re * resultComplex(it).re) / ((h*w) * (h*w));

      //computeEverything : finds the peak, calculates sum, and then returns the PCE
      corrMatrix(i,j) = computeEverything(h, w, result);
    }
    corrTimer.stop();
    // Increment by the number of correlations in 1 batch

    // Stopping the timer because the beginning of this loop is just initializations
    overallTimer.stop();
    idx += num;
    cleanup_threads();
    cleanup();
  }
  
  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer.elapsed(), "s");
  writeln("PRNU Time: ", prnuTimer.elapsed(), "s");
  writeln("FFT Time: ", fftTimer.elapsed(), "s");
  writeln("Cross Time : ", crossTimer.elapsed(), "s");
  writeln("Corr TIme : ", corrTimer.elapsed(), "s");

  flushWriteln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");
  flushWriteln("Sleep before destoy plans for prnu data");
  sleep(60);
  forall i in numDomain {
    prnuDestroy(data(i));
  }
  complete_destruction();
  flushWriteln("Sleep after destroy plans for prnu data");

  /* IMP : Seeing a memory leak of 4GB when run for 10 images with batch size of 10 */
  sleep(60);
  printGlobalMemory("******* After everything is done ******");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}