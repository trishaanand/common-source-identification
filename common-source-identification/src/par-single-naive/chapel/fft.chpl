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

proc computeEverything(h : int, w : int, prnuComplex : [] complex, prnuRotComplex : [] complex) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    var resultComplex : [imageDomain] complex;
    
    // Calculate the point wise product of both matrices
    for (i,j) in imageDomain do {
        resultComplex(i,j) = prnuComplex(i,j) * prnuRotComplex(i,j);
    }

    // Calculate inverse FFT of the result
    calculateFFT(resultComplex, FFTW_BACKWARD);
    
    // Scale the result
    resultComplex /= h*w;
    
    var result : [imageDomain] real;

    var max, sum : real;
    sum = 0.0;
    var maxI, maxJ : int;
    for (i,j) in imageDomain do {
        result(i,j) = resultComplex(i,j).re;
        
        //Find the peak
        if (max < result(i,j)) {
            max = result(i,j);
            maxI = i;
            maxJ = j;
        }
    }
    
    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w) then w else maxJ + 5;

    for (i,j) in imageDomain do {
        if(!(i > lowI && i < highI && j > lowJ && j < highJ )) {
            sum += result(i,j) * result(i,j);
        }
    }
    
    //Calculate average energy
    var energy : real;
    // energy = sum/((h*w) - ((highI-lowI + 1)*(highJ-lowJ + 1)));
    energy = sum/((h*w) - 121);
    var PCE = (max * max) / energy;
    return PCE;
}
