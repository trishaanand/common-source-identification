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
config const num : int = 100; 

var data : prnu_data;  

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnuComplex(h : int, w : int, image : [] RGB) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};

  var prnu : [imageDomain] real;
  var prnuComplex : [imageDomain] complex; 
  
  prnuExecute(prnu, image, data);

  for(i,j) in imageDomain {
    prnuComplex(i,j) = prnu(i,j) + 0i;
  }

  return prnuComplex;
}

proc rotated180Prnu(h : int, w : int, prnu : [] complex) {
  const imageDomain: domain(2) = {0..#h, 0..#w};
  var prnuRot : [imageDomain] complex;

  /* Rotate the matrix 180 degrees */
  for (i,j) in imageDomain do 
    prnuRot(i,j) = prnu(h-i-1, w-j-1);

  return prnuRot;
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

  var overallTimer, fftTimer, corrTimer : Timer;
  var t1Timer, t2Timer, t3Timer, t4Timer, t5Timer : real;
  var sumt1Timer, sumt2Timer, sumt3Timer, sumt4Timer, sumt5Timer : real;

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const numDomain : domain(1) = {1..n};

  flushWriteln("Running Common Source Identification...");
  flushWriteln("  ", n, " images");
  flushWriteln("  ", numLocales, " locale(s)");

  /* ************************* Start here ********************* */
  for (i,j) in corrDomain {
    if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
    }
  }

  /* Perform all initializations here */
  prnuInit(h, w, data);
  
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
    
    for it in localSubDom {
      var (i,j) = crossTuples(it);
      imgSparseDom += i;
      imgSparseDom += j;
      prnuSparseDom += i; 
      prnuRotSparseDom += j;
    }
    var images : [imgSparseDom][imageDomain] RGB;
    var prnuArray : [prnuSparseDom][imageDomain] complex;
    var prnuRotArray : [prnuRotSparseDom][imageDomain] complex;
    overallTimer.stop();

    for i in imgSparseDom {
      readJPG(images[i], imageFileNames[i]);
    }

    /* Start the timer here */
    overallTimer.start();

    for i in imgSparseDom {
      // Calculate the prnu of the required images
      var prnu = calculatePrnuComplex(h, w, images[i]);
      if prnuSparseDom.member(i) {
        prnuArray(i) = prnu;
        calculateFFT(prnuArray(i), FFTW_FORWARD);
      }
      if prnuRotSparseDom.member(i) {
        prnuRotArray(i) = rotated180Prnu(h, w, prnu);
        calculateFFT(prnuRotArray(i), FFTW_FORWARD);
      }
    }

    /* Calculate correlation now */
    for it in localSubDom {
      var (i,j) = crossTuples(it);
      corrMatrix(i,j) = computeEverything(h, w, prnuArray(i), prnuRotArray(j));
    }

    // Increment by the number of correlations in 1 batch
    idx += num;

    // Stopping the timer because the beginning of this loop is just initializations
    overallTimer.stop();
  }
  
  prnuDestroy(data);

  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer.elapsed(), "s");

  flushWriteln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}
