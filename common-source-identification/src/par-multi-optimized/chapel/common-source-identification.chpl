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

use Memory;

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
      if (i > j) {
        channel.writef("%{#####.#####} ", array[i, j]);
      } else if (i < j) {
        channel.writef("%{#####.#####} ", array[j, i]);
      } else {
        channel.writef("%{#####.#####} ", 0.0);
      }
    }
    channel.writeln();
  }
}

proc complexAndRotate(i : int, j : int, h : int, w : int, 
                      prnu : [] real, prnuComplex : [] complex, prnuRotComplex : [] complex) {
  prnuComplex(i,j) = prnu(i,j) + 0i;
  prnuRotComplex(i,j) = prnu(h-i-1, w-j-1) + 0i;
}

proc myWrite(s...?) {
  stdout.writeln(s);
  stdout.flush();
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
  // const totMem = + reduce Locales.physicalMemory();
  // writeln("Total mem available across all locales : ", totMem);

  // forall (i,j) in corrDomain {
  //   if (i > j) {
  //     var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
  //     crossTuples(idx) = (i-1,j-1);
  //   }
  // }

  forall loc in Locales {
    on loc do {
      const x = here.physicalMemory(unit=MemUnits.GB, retType=real);
      writeln("On locale: ", loc.id, " memory: ", x);
      // writeln("On locale: ", loc.id, " memory: ", locale.physicalMemory(MemUnits.MB, x));

      // var subNumDomainPRNU : sparse subdomain(numDomain);
      // var subNumDomainROT : sparse subdomain(numDomain);
      // var arr : [subNumDomainPRNU] real; 
      // for idx in crossDomain.localSubdomain() {
      //   var (x,y)= crossTuples(idx);
      //   subNumDomainPRNU += x;
      //   subNumDomainROT += y;
      //   writeln("subNumDomainPRNU IRV: ", arr.IRV);
      //   writeln("On locale ", here.id, " tuple stored in the sparse domain is : ", crossTuples(idx));
      // }
    }
  }


  writeln("End of execution");
}

proc printGlobalMemory(s : string) {
  coforall loc in Locales do on loc {
    myWrite(s , " on locale: ", here.id, " mem: ", memoryUsed()/1000000);
  }
}

proc printLocalMemory(s : string) {
  myWrite(s , " on locale: ", here.id, " mem: ", memoryUsed()/1000000);
}

/**
  TODO: Potential improvements in the codebase:
  1. Figure out why removing the prnuTemp array & using flags doesn't work on the local machines. 
  2. CorrMatrix should also be block mapped. The println for the CorrMatrix is outside the timer loop.
  */

proc run() {
  printGlobalMemory("Before doing anything");
  
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
  var crossNum = n * (n-1) / 2;
  const tupleCrossDomain = {0..#crossNum} dmapped Block({0..#crossNum}); 

  // Define all the classes here
  class Im {
    var images : [numDomain][imageDomain] RGB;
    var data : [numDomain] prnu_data;  
  }
  
  printGlobalMemory("Before PRNU");
  
  var prnuArray, prnuRotArray : [numDomain][imageDomain] complex;
  printGlobalMemory("After PRNU");
  var overallTimer, prnuTimer, fftTimer, corrTimer, crossTimer, copyTimer : Timer;

  const replCrossDomain = {0..#crossNum} dmapped Replicated();
  var crossTuples : [replCrossDomain] 2*int;
  
  myWrite("Running Common Source Identification...");
  myWrite("  ", n, " images");
  myWrite("  ", numLocales, " locale(s)");
  
  for loc in Locales do on loc {
    forall (i,j) in corrDomain {
      if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
      }
    }
    printLocalMemory("Before reading images");
  }
  
  var imObj = new unmanaged Im ();

  for i in numDomain {
    var image : [imageDomain] RGB;
    readJPG(image, imageFileNames[i]);
    imObj.images[i] = image; //Assumption that this would run on the locale image(i) belongs to.
  }

  /* Perform all the initializations */
  coforall loc in Locales do on loc {
    for i in numDomain.localSubdomain() {
      prnuInit(h,w,imObj.data(i));
    }
    printLocalMemory("After creating class storing image and data");
  }

  /* Start the timer to measure the ops */
  overallTimer.start();
  
  coforall loc in Locales do on loc {
    var prnus : [numDomain.localSubdomain()][imageDomain] real;
    forall i in numDomain.localSubdomain() {
      prnuExecute(prnus[i], imObj.images[i], imObj.data[i]);
      forall (k, j) in imageDomain {
        complexAndRotate(k, j, h, w, prnus[i], prnuArray[i], prnuRotArray[i]);
      }
    }
  }

  delete imObj;
  printGlobalMemory("After deleting class obj");
  
  coforall loc in Locales do on loc {
    //Create plans for the FFT for the prnu arrays : MUST BE SINGLE THREADED
    { 
      var fftPlanNormal, fftPlanRotate : [numDomain.localSubdomain()] fftw_plan;
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
      if(loc.id != prnuArray[i].locale.id) {
        subNumDomainPRNU += i;
      }
      if(loc.id != prnuRotArray[j].locale.id) {
        subNumDomainROT += j;
      }
    }
    printLocalMemory("Before alloc resultComplex");
    var resultComplex : [tupleCrossDomain.localSubdomain()][imageDomain] complex;
    
    {
      var prnuRemote : [subNumDomainPRNU][imageDomain] complex;
      printLocalMemory("After prnuRemote");
      var prnuRotRemote : [subNumDomainROT][imageDomain] complex;
      printLocalMemory("After prnuRotRemote");
      // // To ensure that the prnu & prnuRot data is not copied multiple times
      var flags, flagsRot : [numDomain] bool; 
      forall idx in tupleCrossDomain.localSubdomain() {
        var prnuTemp, prnuRotTemp : [imageDomain] complex;
        var (i,j) = crossTuples(idx);
        if(loc.id != prnuArray[i].locale.id) {
          if (flags[i] == false) {
            prnuRemote[i] = prnuArray[i]; //copy from locale which stores prnuArray[i] to store locally
            flags[i] = true;
          } 
          prnuTemp = prnuRemote[i];
        } else {
          prnuTemp = prnuArray[i];
        }
        if(loc.id != prnuRotArray[j].locale.id) {
          if (flagsRot[j] == false) {
            prnuRotRemote[j] = prnuRotArray[j];
            flagsRot[j] = true;
          } 
          prnuRotTemp = prnuRotRemote[j];
        } else {
          prnuRotTemp = prnuRotArray[j];
        }
        
        resultComplex(idx) = prnuTemp * prnuRotTemp;
      }
    }
    
    printLocalMemory("After rotRemote & remote are deleted");
    {
      var fftPlanBack : [tupleCrossDomain.localSubdomain()] fftw_plan;
      
      // Plan all the FFTs for the resultComplex in a serialized fashion
      for idx in tupleCrossDomain.localSubdomain() {
        fftPlanBack(idx) = planFFT(resultComplex(idx), FFTW_BACKWARD) ; 
      }
      forall idx in tupleCrossDomain.localSubdomain() {
        execute(fftPlanBack(idx));
      }
    }
    
    /* Calculate the enery in the resultComplex array */
    forall idx in tupleCrossDomain.localSubdomain() {
      //Save the real part of result array, scale and square it.
      var result = (resultComplex(idx).re * resultComplex(idx).re) / ((h*w) * (h*w)); 
    
      var (i,j) = crossTuples(idx);
      corrMatrix(i,j) = computeEverything(h, w, result);
    }
  }
  
  
  overallTimer.stop();
  
  printGlobalMemory("After everything is done");
  
  myWrite("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  myWrite("Time: ", overallTimer.elapsed(), "s");

  var nrCorrelations = (n * (n - 1)) / 2;
  myWrite("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    myWrite("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }
  
  myWrite("End");
}