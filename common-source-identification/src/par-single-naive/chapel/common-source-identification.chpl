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
proc calculatePrnuComplex(h : int, w : int, imageFileName : string) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};
  var image : [imageDomain] RGB;

  /* allocate a prnu_data record */
  var data : prnu_data;  
  var prnu : [imageDomain] real;

  /* Read in the first image. */
  readJPG(image, imageFileName);

  prnuInit(h, w, data);
  prnuExecute(prnu, image, data);
  prnuDestroy(data);

  var prnuComplex = [ij in imageDomain] prnu(ij) + 0i;

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

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const prnuDomain : domain(1) = {1..n};

  var prnuArray, prnuRotArray : [prnuDomain][imageDomain] complex;
  
  var overallTimer : Timer;

  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");

  overallTimer.start();

  for i in prnuDomain {
    prnuArray(i) = calculatePrnuComplex(h, w, imageFileNames[i]);
    calculateFFT(prnuArray(i), FFTW_FORWARD);
    prnuRotArray(i) = rotated180Prnu(h, w, prnuArray(i));
    calculateFFT(prnuRotArray(i), FFTW_FORWARD);
  }

  /* Calculate correlation now */
  for (i, j) in corrDomain {
    // Only calculating for values below the diagnol of the matrix. The upper half can simply be equated
    // to the lower half
    if(i > j) {
      //call function here.
      corrMatrix(i,j) = computeEverything(h, w, prnuArray(i), prnuRotArray(j));
      corrMatrix(j,i) = corrMatrix(i,j);
    }        
  }

  overallTimer.stop();

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}