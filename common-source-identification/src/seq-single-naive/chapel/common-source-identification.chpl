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

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnuComplex(h : int, w : int, image : [] RGB, prnuComplex : [] complex, ref data : prnu_data) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};

  var prnu : [imageDomain] real;
  
  prnuExecute(prnu, image, data);

  for(i,j) in imageDomain {
    prnuComplex(i,j) = prnu(i,j) + 0i;
  }
}

proc rotated180Prnu(h : int, w : int, prnu : [] complex, prnuRot : [] complex) {
  const imageDomain: domain(2) = {0..#h, 0..#w};

  /* Rotate the matrix 180 degrees */
  for (i,j) in imageDomain do 
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

  var corrMatrix : [corrDomain] real;
  var nrCorrelations = (n * (n - 1)) / 2;
  var crossDomain : domain(1) = {0..#nrCorrelations};
  var crossTuples : [crossDomain] 2*int;
  
  var data : prnu_data;  
  var prnu, prnuRot, resultComplex : [imageDomain] complex;

  var overallTimer : Timer;


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

  var fwPlan = plan_dft(prnu, prnu, FFTW_FORWARD, FFTW_ESTIMATE);
  var fwPlanRot = plan_dft(prnuRot, prnuRot, FFTW_FORWARD, FFTW_ESTIMATE);
  var bwPlan = plan_dft(resultComplex, resultComplex, FFTW_BACKWARD, FFTW_ESTIMATE);

  for idx in crossDomain {
    var (i,j) = crossTuples(idx);

    // Read both the images from disk 
    var images : [2][imageDomain] RGB;
    readJPG(images[0], imageFileNames[i]);
    readJPG(images[1], imageFileNames[j]);
    
    overallTimer.start();

    var prnuTemp : [imageDomain] complex;
    calculatePrnuComplex(h, w, images[0], prnu, data);
    calculatePrnuComplex(h, w, images[1], prnuTemp, data);
    rotated180Prnu(h, w, prnuTemp, prnuRot);

    // Calculate the FFT on the prnu & rot arrays
    execute(fwPlan);
    execute(fwPlanRot);

    // Calculate the dot product
    for (x,y) in imageDomain {
      resultComplex(x,y) = prnu(x,y) * prnuRot(x,y);
    }

    // Inverse FFT
    execute(bwPlan);

    corrMatrix(i,j) = computePCE(h, w, resultComplex);
    overallTimer.stop();
  }

  // Cleanup everything
  prnuDestroy(data);
  destroy_plan(fwPlan);
  destroy_plan(fwPlanRot);
  destroy_plan(bwPlan);
  cleanup();

  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer.elapsed(), "s");

  flushWriteln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}
