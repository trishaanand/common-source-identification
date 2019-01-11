use FFTW;
use Time;

class PlusReduceOp: ReduceScanOp {
  type eltType;
  var  value: eltType;
  proc identity         return 0: eltType;
  proc accumulate(elm)  { value = value + elm; }
  proc accumulateOntoState(ref state, elm) { state = state + elm; }
  proc combine(other)   { value = value + other.value; }
  proc generate()       return value;
  proc clone()          return new unmanaged PlusReduceOp(eltType=eltType);
}

/* prnu is passed by reference. Hence not returning any value from the FFT calculation */ 
proc calculateFFT(prnu : [] complex, sign : c_int) {
    var fftPlan = plan_dft(prnu, prnu, sign, FFTW_ESTIMATE);
    execute(fftPlan);
}

proc planFFT(prnu : [] complex, sign : c_int) {
    return plan_dft(prnu, prnu, sign, FFTW_ESTIMATE);
}

//Number of results : N*(N-1)/2
proc computeEverything(h : int, w : int, prnuComplex : [] complex, prnuRotComplex : [] complex) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    var resultComplex : [imageDomain] complex;
    
    var t1Timer, t2Timer, t3Timer, t4Timer, t5Timer : Timer;

    // writeln("In computeEverthign, (100,100): ", prnuComplex(100,100));
    // writeln("In computeEverthign rotated, (100,100): ", prnuRotComplex(100,100));
    // Calculate the point wise product of both matrices
    
    t1Timer.start();
    resultComplex = prnuComplex * prnuRotComplex;
    t1Timer.stop();
    
    // writeln("after cross correlation, (100,100): ", resultComplex(100,100));

    // Calculate inverse FFT of the result
    t2Timer.start();
    calculateFFT(resultComplex, FFTW_BACKWARD);
    t2Timer.stop();

    var result : [imageDomain] real;
    var max : real;
    var maxI, maxJ : int;
    
    // Scale the result and square the real component
    t3Timer.start();
    result = resultComplex.re;
    result = (result * result) / ((h*w) * (h*w)); 
    t3Timer.stop();

    t4Timer.start();
    var (maxVal, maxLoc) = maxloc reduce zip(result, imageDomain);
    max = maxVal;
    maxI = maxLoc(1);
    maxJ = maxLoc(2);
    t4Timer.stop();

    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w) then w else maxJ + 5;

    t5Timer.start();
    // var innerDomain : domain(2) = {lowI..highI-1, lowJ..highJ-1};
    // var ignoreSum = + reduce result[innerDomain];
    // // writeln("Ignore sum is ", ignoreSum);
    // var sum = + reduce result;
    // // writeln("Total sum via reduce is ", sum, "and sum after ignoring is : ", (sum - ignoreSum));
    // sum = sum - ignoreSum;
    var sum=0.0;
    // var ignoreSum1 =0.0;
    for (i,j) in imageDomain {
        if(!(i > lowI && i < highI && j > lowJ && j < highJ )) {
            sum += result(i,j);
        // } else {
        //     ignoreSum1 += result(i,j);
        }
    }
    t5Timer.stop();

    //Calculate average energy
    var energy : real;
    // energy = sum/((h*w) - ((highI-lowI + 1)*(highJ-lowJ + 1)));
    energy = sum/((h*w) - 121);
    var PCE = (max) / energy;
    return (PCE, t1Timer.elapsed(), t2Timer.elapsed(), t3Timer.elapsed(), t4Timer.elapsed(), t5Timer.elapsed());
}
