/* Use assertions. */
use Assert;

/* Use sorting. */
use Sort;

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

use VisualDebug;

/* Configuration parameters */
config const imagedir : string = "images";
config const writeOutput : bool = false;
var data : prnu_data;  

/* Add a directory to a file name. */
proc addDirectory(fileName : string, dir : string) : string {
  return dir + "/" + fileName;
}

/* Get an array of the file names of the images in directory dir. */
proc getImageFileNames(dir : string) {
    var imageFiles = listdir(dir);
    sort(imageFiles);
    return addDirectory(imageFiles, dir);
}

/* Write a real array to a file. */
proc write2DRealArray(array : [] real, fileName :string) {
  assert(array.rank == 2);
  var file = open(fileName, iomode.cw);
  var channel = file.writer();

  for i in array.domain.dim(1) {
    for j in array.domain.dim(2) {
      channel.writef("%{#####.#####} ", array[i, j]);
    }
    channel.writeln();
  }
}

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnu(h : int, w : int, image : [] RGB, prnuComplex : [] complex, prnuRotComplex : [] complex) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};

  /* allocate a prnu_data record */
  // var data : prnu_data;  
  var prnu : [imageDomain] real;

  // prnuInit(h, w, data);
  prnuExecute(prnu, image, data);
  // prnuDestroy(data);

  forall (i, j) in imageDomain {
    complexAndRotate(i, j, h, w, prnu, prnuComplex, prnuRotComplex);
  }
}

proc complexAndRotate(i : int, j : int, h : int, w : int, 
                      prnu : [] real, prnuComplex : [] complex, prnuRotComplex : [] complex) {
  prnuComplex(i,j) = prnu(i,j) + 0i;
  prnuRotComplex(i,j) = prnu(h-i-1, w-j-1) + 0i;
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
  var corrMatrix : [corrDomain] real;

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const numDomain : domain(1) = {1..n};
  var images : [numDomain][imageDomain] RGB;
  var t1Timer, t2Timer, t3Timer, t4Timer, t5Timer : real;
  var sumt1Timer, sumt2Timer, sumt3Timer, sumt4Timer, sumt5Timer : real;

  // var data : [numDomain] prnu_data;  
  var prnuArray, prnuRotArray : [numDomain][imageDomain] complex;
  
  var overallTimer, fftTimer, corrTimer : Timer;

  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  writeln("  ", here.numPUs(false, true), " cores");
  writeln("  ", here.maxTaskPar, " maxTaskPar");

  /* Perform all the initializations */
  prnuInit(h,w,data);
  forall i in numDomain {
    readJPG(images[i], imageFileNames[i]);
    
  }

  /* Start the timer to measure the ops */
  overallTimer.start();
  
  fftTimer.start();
  for i in numDomain {
    calculatePrnu(h, w, images[i], prnuArray(i), prnuRotArray(i));
    calculateFFT(prnuArray(i), FFTW_FORWARD);
    calculateFFT(prnuRotArray(i), FFTW_FORWARD);
  }

  writeln("**************** hohohoho *********");

  fftTimer.stop();

  /* Calculate correlation now */
  corrTimer.start();
  for (i, j) in corrDomain {
    // Only calculating for values below the diagnol of the matrix. The upper half can simply be equated
    // to the lower half
    if(i < j) {
      //call function here.
      // writeln("************** i,j: ", i, " ", j);
      (corrMatrix(i,j), t1Timer, t2Timer, t3Timer, t4Timer, t5Timer) = computeEverything(h, w, prnuArray(i), prnuRotArray(j));
      corrMatrix(j,i) = corrMatrix(i,j);
      // sumt1Timer += t1Timer;
      // sumt2Timer += t2Timer;
      // sumt3Timer += t3Timer;
      // sumt4Timer += t4Timer;
      // sumt5Timer += t5Timer;
    }
  }
  corrTimer.stop();
  writeln("**************** hahahahaha *********");
  
  overallTimer.stop();
  // forall i in numDomain {
  //   prnuDestroy(data[i]);
  // }
    prnuDestroy(data);
  

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  writeln("PRNU + FFT Time: ", fftTimer.elapsed(), "s");
  writeln("Corr TIme : ", corrTimer.elapsed(), "s");
  // writeln("T1 Timer: ", sumt1Timer, "s");
  // writeln("T2 Timer: ", sumt2Timer, "s");
  // writeln("T3 Timer: ", sumt3Timer, "s");
  // writeln("T4 Timer: ", sumt4Timer, "s");
  // writeln("T5 Timer: ", sumt5Timer, "s");

  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}