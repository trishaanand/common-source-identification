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

/* Configuration parameters */
config const imagedir : string = "images";
config const writeOutput : bool = false;
var data : prnu_data;  

/* Add a directory to a file name. */
proc addDirectory(fileName : string, dir : string) : string {
  return dir + "/" + fileName;
}

/* Get an array of the file names of the images in directory dir. */
proc getImageFileNamesFullPath(dir : string) {
    var imageFiles = listdir(dir);
    sort(imageFiles);
    return addDirectory(imageFiles, dir);
}

proc getRotatedFilename(fileName: string) {
  return "rot-" + fileName;
}

/* Get an array of the file names of the images in directory dir. */
proc getImageFileNames(dir : string) {
    var imageFiles = listdir(dir);
    sort(imageFiles);
    return imageFiles;
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
proc calculatePrnuComplex(h : int, w : int, image : [] RGB) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};

  /* allocate a prnu_data record */
  
  var prnu : [imageDomain] real;
  var prnuComplex : [imageDomain] complex; 
  
  // prnuInit(h, w, data);
  prnuExecute(prnu, image, data);
  // prnuDestroy(data);

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
  var imageFileNames = getImageFileNamesFullPath(imagedir);

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
  var images : [numDomain][imageDomain] RGB;

  var prnuArray, prnuRotArray : [numDomain][imageDomain] complex;
  
  var overallTimer, fftTimer, corrTimer : Timer;

  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");

  /* Perform all initializations here */
  prnuInit(h, w, data);
  for i in numDomain {
    readJPG(images[i], imageFileNames[i]);
  }

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

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  writeln("PRNU + FFT Time: ", fftTimer.elapsed(), "s");
  writeln("Corr TIme : ", corrTimer.elapsed(), "s");
  writeln("T1 Timer: ", sumt1Timer, "s");
  writeln("T2 Timer: ", sumt2Timer, "s");
  writeln("T3 Timer: ", sumt3Timer, "s");
  writeln("T4 Timer: ", sumt4Timer, "s");
  writeln("T5 Timer: ", sumt5Timer, "s");

  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}
