use FFTW_MT;
use Math;

/* prnu is passed by reference. Hence not returning any value from the FFT calculation */ 
proc calculateFFT(prnu : [] complex, sign : c_int) {
    // Potentially use plan_dft_r2c which converts real arrays to complex arrays during the FFT conversion
    // Check: https://chapel-lang.org/docs/modules/packages/FFTW.html#module-FFTW
    // Currently, we are using complex to complex out-of-place DFT plans.
    var fftPlan = plan_dft(prnu, prnu, sign, FFTW_ESTIMATE);
    execute(fftPlan);
}

proc calculateFFTR2C(prnu : [] real, prnuComplex : [] complex) {
    var fftPlan = plan_dft_r2c(prnu, prnuComplex, FFTW_ESTIMATE);
    execute(fftPlan);
}

proc calculateFFTC2R(prnuComplex : [] complex, prnu : [] real) {
    var fftPlan = plan_dft_c2r(prnuComplex, prnu, FFTW_ESTIMATE);
    execute(fftPlan);
}

proc computeEverything(h : int, w : int, prnuComplex : [] complex, prnuRotComplex : [] complex, fftDomain : domain) {
    const imageDomain: domain(2) = {0..#h, 0..#w};

    var resultComplex : [fftDomain] complex;
    var result : [fftDomain] real;

    for (i,j) in imageDomain {
        resultComplex(i,j,0) = prnuComplex(i,j,0) * prnuRotComplex(i,j,0);
    }
    // resultComplex = prnuComplex * prnuRotComplex;

    calculateFFTC2R(resultComplex, result);
    
    // Scale the result
    result /= (h*w*2);
    
    var max, sum : real;
    sum = 0.0;
    var maxI, maxJ : int;
    for (i,j) in imageDomain do {
        //Find the peak
        if (max < result(i,j,0)) {
            max = result(i,j,0);
            maxI = i;
            maxJ = j;
        }
    }
    
    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w/2) then w else maxJ + 5;

    for (i,j) in imageDomain do {
        if(!(i > lowI && i < highI && j > lowJ && j < highJ )) {
            sum += result(i,j,0) * result(i,j,0);
        }
    }
    
    //Calculate average energy
    var energy : real;
    // energy = sum/((h*w) - ((highI-lowI + 1)*(highJ-lowJ + 1)));
    energy = sum/((h*w) - 121);
    // writeln("energy: ", energy);
    var pce = (max * max) / energy;
    // writeln("pce: ", pce);
    return pce;
}
