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

/* Configuration parameters */
config const imagedir : string = "images";
config const writeOutput : bool = false;

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

  var overallTimer : Timer;

  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");

  overallTimer.start();

  /* Below is example code that contains part of the pieces that you need to 
   * compute the correlation matrix.
   * 
   * It shows how to read in a JPG image, how to compute a noise pattern.
   * Modify the code to compute the correlation matrix.
   */
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h,0..#w};
  var image : [imageDomain] RGB;

  /* Read in the first image. */
  readJPG(image, imageFileNames.front());

  /* allocate a prnu_data record */
  var data : prnu_data;
  var prnu : [imageDomain] real;

  prnuInit(h, w, data);
  prnuExecute(prnu, image, data);
  prnuDestroy(data);

  overallTimer.stop();

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
    // for now, also write the prnu noise pattern, can be removed
    write2DRealArray(prnu, "prnu");
  }
}