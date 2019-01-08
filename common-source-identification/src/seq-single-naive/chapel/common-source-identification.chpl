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
proc calculatePrnu(h : int, w : int, imageFileName : string) {
  
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

  return prnu;
}

proc rotated180Prnu(h : int, w : int, prnu : [] real) {
  const imageDomain: domain(2) = {0..#h, 0..#w};
  var prnuRot : [imageDomain] real;

  /* Rotate the matrix 180 degrees */
  for (i,j) in imageDomain do 
    prnuRot(i,j) = prnu(h-i-1, w-j-1);

  return prnuRot;
}

proc main() {
  /* Obtain the images. */
  // images/Pentax_XXX.JPG
  writeln("Start here");
  var imageFileNames = getImageFileNamesFullPath(imagedir);
  // Arpit - Getting the list of files to write prnu for an image with the same filename as the original image
  // Can also be derived by splitting the complete image but this was just easier.
  var imageFilesPRNU = getImageFileNames(imagedir);

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
  
  var prnuArray, prnuRotArray : [prnuDomain][imageDomain] real;

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

  // Calculating the prnus and really trying hard to save to an array
  var prnuI, prnuRotJ : [imageDomain] real;

  for i in prnuDomain {
    var prnu = calculatePrnu(h, w, imageFileNames[i]);
    var prnuRot = rotated180Prnu(h, w, prnu);
    if(i == 1) {  
      prnuI = prnu;

      writeln("Before fft prnuOrig (100,100) : ", prnu(100,100));
      writeln("Before fft prnuI (100,100) : ", prnuI(100,100));
    } else if (i == 2) {
      prnuRotJ = prnuRot;
      writeln("Before fft prnuRotOrig (100,100) : ", prnuRot(100,100));
      writeln("Before fft prnuRotJ (100,100) : ", prnuRotJ(100,100));
    }
    // writeln("Before fft (100,100) : ", prnu(100,100));

    // prnuArray[i][.., ..] = prnu;
    // writeln("Before fft for array (100,100) : ", prnuArray[i]);
    // // writeln("Value of i: ", i);
    // var prn : [imageDomain] real;
    // prn = prnuArray[i];
    // writeln("Before fft array (100,100) : ", prn(100, 100));
    // prnuRotArray(i) = prnuRot;

    // if(writeOutput) {
    //   write2DRealArray(prnu, imageFilesPRNU[i]);
    //   write2DRealArray(prnuRot, getRotatedFilename(imageFilesPRNU[i]));
    // }
  }

  computeEverything(h, w, prnuI, prnuRotJ);
  /* Calculate correlation now */
  // for i in prnuDomain {
  //   for j in prnuDomain {
  //     if (i != j) {
  //       //call function here.
  //       // writeln("i,j: ", i, j);
  //       // writeln("Before computeEverything (100,100) : ", prnuArray[i](100, 100));
  //       var prnu = prnuArray(i);
  //       corrMatrix(i,j) = computeEverything(h, w, prnuArray(i), prnuRotArray(j));
  //     }
  //   }
  // }

  overallTimer.stop();

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
    // for now, also write the prnu noise pattern, can be removed
    // write2DRealArray(prnu, imageFiles.front());
  }

  writeln("End");
}