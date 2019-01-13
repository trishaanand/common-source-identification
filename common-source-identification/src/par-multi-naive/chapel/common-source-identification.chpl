/* Use assertions. */
use Assert;

/* Use sorting. */
use Sort;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

use BlockDist;
use ReplicatedDist;

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
proc calculatePrnu(h : int, w : int, image : [] RGB, prnu : [] real, prnuComplex : [] complex, prnuRotComplex : [] complex, ref data : prnu_data) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};

  prnuExecute(prnu, image, data);

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
  // tryRun();
  run();
}


proc tryRun() {
  /* Obtain the images. */
  var imageFileNames = getImageFileNames(imagedir);

  /* n represents the number of images that have to be correlated. */
  var n = imageFileNames.size;

  /* h, w will represent the height and width of an image or PRNU noise pattern 
   * throughout the code.
   */
  var h, w : int;
  (h, w) = getDimensionsJPG(imageFileNames.front());

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const numDomain = {1..n} dmapped Block({1..n});
  
  var images : [numDomain][imageDomain] RGB;
  
  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  
  // startVdebug("vdata");
  for i in numDomain {
    var image : [imageDomain] RGB;
    on Locales[i % numLocales] do {
      // writeln("i : ", i, " on locale : ", here.id);
    }
  }
  // stopVdebug();
  // writeln("End of execution");
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
  const numDomain = {1..n} dmapped Replicated();
  var crossNum = n * (n-1) / 2;
  const crossDomain = {0..#crossNum} dmapped Replicated();
  
  var images : [numDomain][imageDomain] RGB;
  var fftPlanNormal, fftPlanRotate : [numDomain] fftw_plan;
  var fftPlanBack : [crossDomain] fftw_plan;

  var t1Timer, t2Timer, t3Timer, t4Timer, t5Timer : real;
  var sumt1Timer, sumt2Timer, sumt3Timer, sumt4Timer, sumt5Timer : real;

  var prnuArray, prnuRotArray : [numDomain][imageDomain] complex;
  var data : [numDomain] prnu_data;  
  var prnus : [numDomain][imageDomain] real;
  var overallTimer, prnuTimer, fftTimer, corrTimer, crossTimer : Timer;

  var resultComplex : [crossDomain][imageDomain] complex;
  var crossTuples : [crossDomain] 2*int;
  
  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  
  
  for i in numDomain {
    var image : [imageDomain] RGB;
    readJPG(image, imageFileNames[i]);
    coforall loc in Locales do on loc {
      images[i] = image;
    }
  }

  // writeln("After reading all the images on locale 0");
  /* Perform all the initializations */
  coforall loc in Locales do on loc {
    for i in numDomain {
      prnuInit(h,w,data(i));
      // writeln(i, " on prnu location ", prnus(i).locale.id, " here.id is ", here.id," and i.locale id is ", i.locale.id);
    }
  }
  // writeln("Post reading images and prnuInit");

  /* Start the timer to measure the ops */
  overallTimer.start();
  
  prnuTimer.start();
  forall i in numDomain {
    coforall loc in Locales do on loc {
      // writeln("For i ", i, "  On locale ", here.id);
      prnuExecute(prnus[i], images[i], data[i]);
      // writeln("For i ", i, " on locale ", here.id, " post prnuExecute ");
      forall (k, j) in imageDomain {
        complexAndRotate(k, j, h, w, prnus[i], prnuArray[i], prnuRotArray[i]);
      }
    }
    // calculatePrnu(h, w, images[i], prnus[i], prnuArray(i), prnuRotArray(i), data(i));
  }
  prnuTimer.stop();
  // writeln("Post calculating all the PRNU");
  
  // startVdebug("vdata");
  
  //Create plans for the FFT for the prnu arrays : MUST BE SINGLE THREADED
  fftTimer.start();
  for i in numDomain {
    coforall loc in Locales do on loc {
      // writeln("In FFT Plan For i: ", i, " locale.id" , here.id);
      fftPlanNormal(i) = planFFT(prnuArray(i), FFTW_FORWARD) ;
      fftPlanRotate(i) = planFFT(prnuRotArray(i), FFTW_FORWARD) ;  
    }
  }

  // writeln("After fft plans for prnu and rotation arrays");
  sync {
    coforall loc in Locales do on loc {
      for i in numDomain {
        // writeln("In FFT execute For i: ", i, " locale.id" , here.id);
        begin {execute(fftPlanNormal(i));}
        begin {execute(fftPlanRotate(i));}
      }
    }
  }
  fftTimer.stop();
  // writeln("After all the fft plans have been executed");
  // stopVdebug();

  // writeln("PRNU Time: ", prnuTimer.elapsed(), "s");
  // writeln("FFT Time: ", fftTimer.elapsed(), "s");

  // Calculate the point wise product of both matrices
  crossTimer.start();
  forall (i, j) in corrDomain {

    var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
    on Locales[idx % numLocales] do {
    // coforall loc in Locales do on loc {
      // Only calculating for values below the diagnol of the matrix. The upper half can simply be equated
      // to the lower half
      if(i > j) {
        //call function here.
        crossTuples(idx) = (i,j);        
        resultComplex(idx) = prnuArray(i) * prnuRotArray(j);
      } 
    }
  }
  crossTimer.stop();
  

  // Plan all the FFTs for the resultComplex in a serialized fashion
  corrTimer.start();
  for idx in crossDomain {
    // coforall loc in Locales do on loc {
    on Locales[idx % numLocales] do {
      fftPlanBack(idx) = planFFT(resultComplex(idx), FFTW_BACKWARD) ; 
    }
  }

  forall idx in crossDomain {
    // coforall loc in Locales do on loc {
    on Locales[idx % numLocales] do {
      execute(fftPlanBack(idx));
    }
  }

  // coforall loc in Locales do on loc {
  forall (idx) in crossDomain {
    on Locales[idx % numLocales] do {
      //Save the real part of result array, scale and square it.
      var result : [imageDomain] real;
      result = resultComplex(idx).re;
      result = (result * result) / ((h*w) * (h*w)); 
    
      //call function here.
      var (i,j) = crossTuples(idx);
      corrMatrix(i,j) = computeEverything(h, w, result);
      corrMatrix(j,i) = corrMatrix(i,j);
    }
  }
  corrTimer.stop();
  
  overallTimer.stop();
  for i in numDomain {
    prnuDestroy(data[i]);
  }
  
  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  writeln("PRNU Time: ", prnuTimer.elapsed(), "s");
  writeln("FFT Time: ", fftTimer.elapsed(), "s");
  writeln("Cross Time : ", crossTimer.elapsed(), "s");
  writeln("Corr TIme : ", corrTimer.elapsed(), "s");

  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}