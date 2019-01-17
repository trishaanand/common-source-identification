/* Use assertions. */
use Assert;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

/* Block distribution module */
use BlockDist;

/* Build the wall! */
use Barriers;

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
config const batchNum : int = 10; 

/* Given a file name this function calculates & returns the prnu data for that image */
proc calculatePrnuComplex(h : int, w : int, image : [] RGB, ref data : prnu_data) {
  
  /* Create a domain for an image and allocate the image itself */
  const imageDomain: domain(2) = {0..#h, 0..#w};
  
  var prnu : [imageDomain] real;
  var prnuComplex : [imageDomain] complex; 
  
  prnuExecute(prnu, image, data);

  forall (i,j) in imageDomain {
    prnuComplex(i,j) = prnu(i,j) + 0i;
  }

  return prnuComplex;
}

proc rotated180Prnu(h : int, w : int, prnu : [] complex) {
  const imageDomain: domain(2) = {0..#h, 0..#w};
  var prnuRot : [imageDomain] complex;

  /* Rotate the matrix 180 degrees */
  forall (i,j) in imageDomain do 
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
  var nrCorrelations = (n * (n - 1)) / 2;
  var crossDomain : domain(1) = {0..#nrCorrelations};

  var crossTuples : [crossDomain] 2*int;
  var totalTime : real;
  var overallTimer : Timer;

  const imageDomain : domain(2) = {0..#h, 0..#w};
  const numDomain : domain(1) = {1..n};
  
  flushWriteln("Running Common Source Identification...");
  flushWriteln("  ", n, " images");
  flushWriteln("  ", numLocales, " locale(s)");
  flushWriteln("  ", batchNum, " batchNum");

  /* ************************* Start here ********************* */
  forall (i,j) in corrDomain {
    if (i > j) {
        var idx = (((i-1) * (i - 1 -1)) / 2 ) + (j - 1);
        crossTuples(idx) = (i,j);
    }
  }

  var idx : int = nrCorrelations - 1;
  
  //batchNum is a config variable which defines max number of correlations per locale in a single batch 
  var num = batchNum * Locales.size;

  while(idx > 0) {
    flushWriteln("");
    printLocalMemory("******* Beginning of while loop ******");
    
    var localRange : int;
    if (idx - num) < 0 {
      localRange = 0;
    } else {
      localRange = idx - num + 1;
    }

    var imgGlobalDom : sparse subdomain(numDomain);
    const localSubDom = {localRange..idx} dmapped Block({localRange..idx});

    for it in localSubDom {
      var (i,j) = crossTuples(it);
      imgGlobalDom += i;
      imgGlobalDom += j;
    }

    var globalImages : [imgGlobalDom][imageDomain]RGB; 
    for i in imgGlobalDom {
      readJPG(globalImages[i], imageFileNames[i]);
    }

    // Define the barrier to delete globalImgs once they've been copied to respective locales
    var b = new Barrier(Locales.size);

    overallTimer.start();

    // Executing the computation on all the locales here
    coforall loc in Locales do on loc {
      printLocalMemory("Beginning of coforall loc ");
      var imgSparseDom, prnuSparseDom, prnuRotSparseDom : sparse subdomain(imgGlobalDom);
      flushWriteln("On locale: ", here.id, ", Local subdomain of correlation matrix is ", localSubDom.localSubdomain());
      for it in localSubDom.localSubdomain() {
        var (i,j) = crossTuples(it);
        imgSparseDom += i;
        imgSparseDom += j;
        prnuSparseDom += i; 
        prnuRotSparseDom += j;
      }
      
      flushWriteln("ImageSparseDom: ", imgSparseDom.size, " on locale: ", here.id, "  ", imgSparseDom);
      flushWriteln("prnuSparseDom: ", prnuSparseDom.size, " on locale: ", here.id, "  ", prnuSparseDom);
      flushWriteln("prnuRotSparseDom: ", prnuRotSparseDom.size, " on locale: ", here.id, "  ", prnuRotSparseDom);
      flushWriteln("On locale: ", here.id, ", the number of correlations for which result is calc is : ", localSubDom.localSubdomain().size);

      printLocalMemory("Before data");
      var data : [numDomain] prnu_data;
      printLocalMemory("After data, before image");  
      var images : [imgSparseDom][imageDomain] RGB;
      printLocalMemory("After images, before PRNU");
      var prnuArray : [prnuSparseDom][imageDomain] complex;
      printLocalMemory("After PRNU, before Rot");
      var prnuRotArray : [prnuRotSparseDom][imageDomain] complex;
      printLocalMemory("After ROT, before fft plan");
      var fftPlan : [prnuSparseDom] fftw_plan;
      printLocalMemory("After fft plan, before fft rot plan");
      var fftRotPlan : [prnuRotSparseDom] fftw_plan;
      printLocalMemory("After fft rot plan, before resultComplex");
      var resultComplex : [localSubDom.localSubdomain()][imageDomain] complex;
      printLocalMemory("After allocating resultComplex");

      // Must be sequential because prnuInit includes an FFT step that cannot be parralellized
      for i in imgSparseDom {
        prnuInit(h, w, data(i));
      }
      forall i in imgSparseDom {
        images[i] = globalImages[i];
      }

      // Call the barrier. Delete the imgGlobalDom
      b.barrier();
      imgGlobalDom.clear();
      flushWriteln("Global Images array size: ", globalImages.size);
      
      flushWriteln("On locale: ", here.id, ", size of images is : ", images.size);

      forall i in imgSparseDom {
        // Calculate the prnu of the required images
        var prnu = calculatePrnuComplex(h, w, images[i], data(i));
        if prnuSparseDom.member(i) {
          prnuArray(i) = prnu;
        }
        if prnuRotSparseDom.member(i) {
          prnuRotArray(i) = rotated180Prnu(h, w, prnu);
        }
      }
      for i in prnuSparseDom {
        fftPlan(i) = planFFT(prnuArray(i), FFTW_FORWARD);  
      }
      for j in prnuRotSparseDom {
        fftRotPlan(j) = planFFT(prnuRotArray(j), FFTW_FORWARD);
      }

      sync {
        begin {
          forall i in prnuSparseDom {
            execute(fftPlan(i));
          }
        }
        begin {
          forall i in prnuRotSparseDom {
            execute(fftRotPlan(i));
          }
        }
      }
      

      forall it in localSubDom.localSubdomain() {
        var (i,j) = crossTuples(it);
        forall (x, y) in imageDomain {
          resultComplex(it)(x,y) = prnuArray(i)(x,y) * prnuRotArray(j)(x,y);
        }
      } 
      var fftPlanBack : [localSubDom.localSubdomain()] fftw_plan;
      
      // Plan the inverse fft
      for it in localSubDom.localSubdomain() {
        fftPlanBack(it) = planFFT(resultComplex(it), FFTW_BACKWARD) ; 
      }

      forall it in localSubDom.localSubdomain() {
        execute(fftPlanBack(it));
      }
      printLocalMemory("Before computing correlation matrix");
      /* Calculate correlation now */
      forall it in localSubDom.localSubdomain() {
        var (i,j) = crossTuples(it);
        
        var result : [imageDomain] real = (resultComplex(it).re * resultComplex(it).re) / ((h*w) * (h*w));

        corrMatrix(i,j) = computeEverything(h, w, result);
      }

      forall i in imgSparseDom {
        prnuDestroy(data(i));
      }
    }
    
    overallTimer.stop();
    // Increment by the number of correlations in 1 batch
    idx -= num;
    printGlobalMemory("At the end of batch processing step in while : ");
  }
  
  flushWriteln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  flushWriteln("Time: ", overallTimer.elapsed(), "s");
  flushWriteln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    flushWriteln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
  }

  flushWriteln("End");
}