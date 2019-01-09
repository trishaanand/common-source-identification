use FFTW_MT;
use Math;

proc computeEverything(h : int, w : int, prnu : [] real, prnuRot : [] real) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    var prnuComplex, prnuRotComplex, resultComplex : [imageDomain] complex;

    writeln("Before fft in fft (100,100) : ", prnu(100,100));
    writeln("Before fft rotated in fft (100,100) : ", prnuRot(100,100));

    for (i,j) in imageDomain do {
        prnuComplex(i,j) = prnu(i,j) + 0i;
        prnuRotComplex(i,j) = prnuRot(i,j) + 0i;
    }

    // TODO: These plans are expensive. Save it for them to be used later
    var forwardPrnu = plan_dft(prnuComplex, prnuComplex, FFTW_FORWARD, FFTW_ESTIMATE);
    var forwardPrnuRot = plan_dft(prnuRotComplex, prnuRotComplex, FFTW_FORWARD, FFTW_ESTIMATE);
    
    execute(forwardPrnu);
    execute(forwardPrnuRot);

    writeln("After fft in fft (100,100) : ", prnuComplex(100,100));
    writeln("After fft rotated in fft (100,100) : ", prnuRotComplex(100,100));

    // Calculate the point wise product of both matrices
    for (i,j) in imageDomain do {
        resultComplex(i,j) = prnuComplex(i,j) * prnuRotComplex(i,j);
    }

    writeln("After multiplication (100,100) : ", resultComplex(100,100));

    // Calculate inverse FFT of the result
    var inverseFFT = plan_dft(resultComplex, resultComplex, FFTW_BACKWARD, FFTW_ESTIMATE);
    execute(inverseFFT);
    
    // Scale the result
    resultComplex /= h*w;
    writeln("after IFFT and scaling (100,100) : ", resultComplex(100,100));
    
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

    writeln("max i: ", maxI, " max j: ", maxJ, " peak: ", max);

    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w) then w else maxJ + 5;

    writeln("lowI: ", lowI, " highI: ", highI, " lowJ: ", lowJ, " highJ: ", highJ);
    for (i,j) in imageDomain do {
        if(!(i > lowI && i < highI && j > lowJ && j < highJ )) {
            sum += result(i,j) * result(i,j);
        }
    }
    
    //Calculate average energy
    var energy : real;
    // energy = sum/((h*w) - ((highI-lowI + 1)*(highJ-lowJ + 1)));
    energy = sum/((h*w) - 121);
    writeln("energy: ", energy);
    var PCE = (max * max) / energy;
    writeln("pce: ", PCE);
    return PCE;
}
