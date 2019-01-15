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
  var t1Timer, t2Timer, t3Timer, t4Timer, t5Timer : real;
  var sumt1Timer, sumt2Timer, sumt3Timer, sumt4Timer, sumt5Timer : real;

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const numDomain : domain(1) = {1..n};

  flushWriteln("Running Common Source Identification...");
  flushWriteln("  ", n, " images");
  flushWriteln("  ", numLocales, " locale(s)");

  /* ************************* Start here ********************* */
  var images : [numDomain][imageDomain] RGB;

  var prnuArray, prnuRotArray : [numDomain][imageDomain] complex;
  
  var overallTimer, fftTimer, corrTimer : Timer;

  /* Perform all initializations here */
  prnuInit(h, w, data);
  for i in numDomain {
    readJPG(images[i], imageFileNames[i]);
  }

  flushWriteln("Read all images");

  overallTimer.start();
  fftTimer.start();

  for i in numDomain {
    prnuArray(i) = calculatePrnuComplex(h, w, images[i]);
    prnuRotArray(i) = rotated180Prnu(h, w, prnuArray(i));
    calculateFFT(prnuArray(i), FFTW_FORWARD);
    calculateFFT(prnuRotArray(i), FFTW_FORWARD);
  }
  fftTimer.stop();

  /* Calculate correlation now */
  corrTimer.start();

  for (i, j) in corrDomain {
    // Only calculating for values below the diagnol of the matrix. The upper half can simply be equated
    // to the lower half
    if(i < j) {
      //call function here.
      (corrMatrix(i,j), t1Timer, t2Timer, t3Timer, t4Timer, t5Timer) = computeEverything(h, w, prnuArray(i), prnuRotArray(j));
      corrMatrix(j,i) = corrMatrix(i,j);
      sumt1Timer += t1Timer;
      sumt2Timer += t2Timer;
      sumt3Timer += t3Timer;
      sumt4Timer += t4Timer;
      sumt5Timer += t5Timer;
    }        
  }

  corrTimer.stop();
  overallTimer.stop();
  prnuDestroy(data);

  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer.elapsed(), "s");
  flushWriteln("PRNU + FFT Time: ", fftTimer.elapsed(), "s");
  flushWriteln("Corr TIme : ", corrTimer.elapsed(), "s");
  flushWriteln("T1 Timer: ", sumt1Timer, "s");
  flushWriteln("T2 Timer: ", sumt2Timer, "s");
  flushWriteln("T3 Timer: ", sumt3Timer, "s");
  flushWriteln("T4 Timer: ", sumt4Timer, "s");
  flushWriteln("T5 Timer: ", sumt5Timer, "s");

  var nrCorrelations = (n * (n - 1)) / 2;
  flushWriteln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}
