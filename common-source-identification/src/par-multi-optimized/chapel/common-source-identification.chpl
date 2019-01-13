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

proc complexAndRotate(i : int, j : int, h : int, w : int, 
                      prnu : [] real, prnuComplex : [] complex, prnuRotComplex : [] complex) {
  prnuComplex(i,j) = prnu(i,j) + 0i;
  prnuRotComplex(i,j) = prnu(h-i-1, w-j-1) + 0i;
}

proc main() {
  // tryRun();
  run();
}

proc printLocales(a) {
  forall i in a do i = i.locale.id;
  writeln(a);
  writeln();
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
  const numDomain : domain(1) = {1..n};
  
  var crossNum = n * (n-1) / 2;
  const crossDomain = {0..#crossNum} dmapped Block({0..#crossNum});
  var testArr : [crossDomain] int;

  printLocales(testArr);
  const corrDomain : domain(2) = {1..n, 1..n};
  var crossTuples : [crossDomain] 2*int;
  
  // var images : [numDomain][imageDomain] RGB;
  
  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  
  forall (i,j) in corrDomain {
    if (i > j) {
      var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
      crossTuples(idx) = (i-1,j-1);
    }
  }

  forall loc in Locales {
    on loc do {
      var subNumDomainPRNU : sparse subdomain(numDomain);
      var subNumDomainROT : sparse subdomain(numDomain);
      var arr : [subNumDomainPRNU] real; 
      for idx in crossDomain.localSubdomain() {
        var (x,y)= crossTuples(idx);
        subNumDomainPRNU += x;
        subNumDomainROT += y;
        writeln("subNumDomainPRNU IRV: ", arr.IRV);
        writeln("On locale ", here.id, " tuple stored in the sparse domain is : ", crossTuples(idx));
      }
    }
  }


  writeln("End of execution");
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
  const numDomain = {1..n} dmapped Block({1..n});
  // const replNumDomain = {1..n} dmapped Replicated();
  var crossNum = n * (n-1) / 2;
  const tupleCrossDomain = {0..#crossNum} dmapped Block({0..#crossNum}); 
  // const crossDomain = {0..#crossNum} dmapped Replicated();
  
  var images : [numDomain][imageDomain] RGB;
  var fftPlanNormal, fftPlanRotate : [numDomain] fftw_plan;
  var fftPlanBack : [tupleCrossDomain] fftw_plan;

  var prnuArray, prnuRotArray : [numDomain][imageDomain] complex;
  var data : [numDomain] prnu_data;  
  var prnus : [numDomain][imageDomain] real;
  var overallTimer, prnuTimer, fftTimer, corrTimer, crossTimer, copyTimer : Timer;

  var resultComplex : [tupleCrossDomain][imageDomain] complex;
  const replCrossDomain = {0..#crossNum} dmapped Replicated();
  var crossTuples : [replCrossDomain] 2*int;
  
  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");
  
  coforall loc in Locales do on loc {
    forall (i,j) in corrDomain {
      if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
      }
    }
  }
  
  for i in numDomain {
    var image : [imageDomain] RGB;
    readJPG(image, imageFileNames[i]);
    images[i] = image; //Assumption that this would run on the locale image(i) belongs to.
  }

  /* Perform all the initializations */
  coforall loc in Locales do on loc {
    for i in numDomain.localSubdomain() {
      prnuInit(h,w,data(i));
    }
  }

  /* Start the timer to measure the ops */
  overallTimer.start();
  
  coforall loc in Locales do on loc {

    forall i in numDomain.localSubdomain() {
      prnuExecute(prnus[i], images[i], data[i]);
      forall (k, j) in imageDomain {
        complexAndRotate(k, j, h, w, prnus[i], prnuArray[i], prnuRotArray[i]);
      }
    }
    //Create plans for the FFT for the prnu arrays : MUST BE SINGLE THREADED
    for i in numDomain.localSubdomain() {
      fftPlanNormal(i) = planFFT(prnuArray(i), FFTW_FORWARD) ;
      fftPlanRotate(i) = planFFT(prnuRotArray(i), FFTW_FORWARD) ;  
    }
    sync {
      forall i in numDomain.localSubdomain() {
        begin {execute(fftPlanNormal(i));}
        begin {execute(fftPlanRotate(i));}
      }
    }
    /*
      1. Copy all the required prnu & prnuRot data to the local machines as required
      2. Calculate the point wise product of both matrices
    */
    /* For reference: 
      const tupleCrossDomain = {0..#crossNum} dmapped Block({0..#crossNum}); 
      crossTuples(idx) = (i, j);
    */
    //Create sparse subdomain to save prnu, and prnuRot arrays from other locales locally.
    var subNumDomainPRNU : sparse subdomain(numDomain);
    var subNumDomainROT : sparse subdomain(numDomain);
    for idx in tupleCrossDomain.localSubdomain() {
      var (i,j)= crossTuples(idx);
      if(loc.id != images[i].locale.id) {
        subNumDomainPRNU += i;
      }
      if(loc.id != images[j].locale.id) {
        subNumDomainROT += j;
      }
    }

    var prnuRemote : [subNumDomainPRNU][imageDomain] complex;
    var prnuRotRemote : [subNumDomainROT][imageDomain] complex;

    // To ensure that the prnu & prnuRot data is not copied multiple times
    var flags, flagsRot : [numDomain] bool; 

// startVdebug("vdata");
    forall idx in tupleCrossDomain.localSubdomain() {
      var prnuTemp, prnuRotTemp : [imageDomain] complex;
      var flagI, flagJ : bool; //To figure out if we should use sparse array, or prnu array for calc result
                                // if true, use prnuArray, if false, use prnuRemote
      var (i,j) = crossTuples(idx);
      if(loc.id != images[i].locale.id) {
        if (flags[i] == false) {
          prnuRemote[i] = prnuArray[i]; //copy from locale which stores prnuArray[i] to store locally
          flags[i] = true;
        } 
        flagI = false;
        prnuTemp = prnuRemote[i];
      } else {
        flagI = true;
        prnuTemp = prnuArray[i];
      }
      if(loc.id != images[j].locale.id) {
        if (flagsRot[j] == false) {
          prnuRotRemote[j] = prnuRotArray[j];
          flagsRot[j] = true;
        } 
        flagJ = false;
        prnuRotTemp = prnuRotRemote[j];
      } else {
        flagJ = true;
        prnuRotTemp = prnuRotArray[j];
      }
      
      // if (flagI && flagJ) {//TT
      //   resultComplex(idx) = prnuArray[i] * prnuRotArray[j];
      // } else if (!flagI && !flagJ) {//FF
      //   resultComplex(idx) = prnuRemote[i] * prnuRotRemote[j];
      // } else if (flagI && !flagJ) {//TF
      //   resultComplex(idx) = prnuArray[i] * prnuRotRemote[j];
      // } else if (!flagI && flagJ) {//FT
      //   resultComplex(idx) = prnuRemote[i] * prnuRotArray[j];
      // }

      resultComplex(idx) = prnuTemp * prnuRotTemp;
    }
    // stopVdebug();
  // Plan all the FFTs for the resultComplex in a serialized fashion
    for idx in tupleCrossDomain.localSubdomain() {
      fftPlanBack(idx) = planFFT(resultComplex(idx), FFTW_BACKWARD) ; 
    }
    forall idx in tupleCrossDomain.localSubdomain() {
      execute(fftPlanBack(idx));
    }
  /* Calculate the enery in the resultComplex array */
    forall idx in tupleCrossDomain.localSubdomain() {
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
    
  
  overallTimer.stop();
  for i in numDomain {
    prnuDestroy(data[i]);
  }
  
  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  // writeln("PRNU Time: ", prnuTimer.elapsed(), "s");
  // writeln("Cross Time : ", crossTimer.elapsed(), "s");
  // writeln("Corr TIme : ", corrTimer.elapsed(), "s");

  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  writeln("End");
}