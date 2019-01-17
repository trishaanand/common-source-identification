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
config const num : int = 50; 

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
  forall (i,j) in corrDomain {
    if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
    }
  }

  var idx : int = 0;

  while(idx < nrCorrelations) {
    overallTimer.start();
    var localRange : int;
    if (idx + num) >= nrCorrelations {
      localRange = nrCorrelations-1;
    } else {
      localRange = idx + num - 1;
    }

    const localSubDom : domain(1) = {idx..localRange};
    var imgSparseDom, prnuSparseDom, prnuRotSparseDom : sparse subdomain(numDomain);
    
    // Sequential iteration because we can't add indices to sparse in forall
    for it in localSubDom {
      var (i,j) = crossTuples(it);
      imgSparseDom += i;
      imgSparseDom += j;
      prnuSparseDom += i; 
      prnuRotSparseDom += j;
    }

    /* Perform all initializations here. Hence stopping the timer */
    overallTimer.stop();
    var data : [numDomain] prnu_data;  
    var images : [imgSparseDom][imageDomain] RGB;
    var prnuArray : [prnuSparseDom][imageDomain] complex;
    var prnuRotArray : [prnuRotSparseDom][imageDomain] complex;
    var fftPlan : [prnuSparseDom] fftw_plan;
    var fftRotPlan : [prnuRotSparseDom] fftw_plan;

    for i in imgSparseDom {
      prnuInit(h, w, data(i));
    }
    for i in imgSparseDom {
      readJPG(images[i], imageFileNames[i]);
    }

    /* Start the timer here */
    overallTimer.start();

    prnuTimer.start();
    forall i in imgSparseDom {
      // Calculate the prnu of the required images
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

    crossTimer.start();
    var resultComplex : [localSubDom][imageDomain] complex;

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

    /* Calculate correlation now */
    forall it in localSubDom {
      var (i,j) = crossTuples(it);
      
      var result : [imageDomain] real = (resultComplex(it).re * resultComplex(it).re) / ((h*w) * (h*w));

      corrMatrix(i,j) = computeEverything(h, w, result);
    }
    corrTimer.stop();
    // Increment by the number of correlations in 1 batch

    // Stopping the timer because the beginning of this loop is just initializations
    overallTimer.stop();
    idx += num;
    forall i in imgSparseDom {
      prnuDestroy(data(i));
    }
  }
  
  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer.elapsed(), "s");
  writeln("PRNU Time: ", prnuTimer.elapsed(), "s");
  writeln("FFT Time: ", fftTimer.elapsed(), "s");
  writeln("Cross Time : ", crossTimer.elapsed(), "s");
  writeln("Corr TIme : ", corrTimer.elapsed(), "s");

  flushWriteln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}