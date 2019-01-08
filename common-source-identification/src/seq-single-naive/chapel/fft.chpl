use FFTW_MT;
use Math;

proc computeEverything(h : int, w : int, prnu : [] real, prnuRot : [] real) {
    const imageDomain: domain(2) = {0..#h,0..#w};
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
    resultComplex /= h*w;
    writeln("after IFFT and scaling (100,100) : ", resultComplex(100,100));
    // TODO: Create new array from the real values in resultComplex
    var result : [imageDomain] real;

    var max, sum : real;
    max = -1.0;
    sum = 0.0;
    var maxI, maxJ : int;
    for (i,j) in imageDomain do {
        result(i,j) = resultComplex(i,j).re;
        
        // Scale the result
        // result(i,j) /= h*w;
        result(i,j) = result(i,j) * result(i,j);
        sum += result(i,j);
        //Find the peak
        if (max < result(i,j)) {
            max = result(i,j);
            maxI = i;
            maxJ = j;
        }

    }
    writeln("max i: ", maxI, " max j: ", maxJ, " peak: ", sqrt(max) );

    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    if ((maxI-5) < 0) {
        lowI = 0;
    } else {
        lowI = maxI - 5;
    }
    if ((maxI + 5) > h) {
        highI = h;
    } else {
        highI = maxI + 5;
    }
    
    if ((maxJ-5) < 0) {
        lowJ = 0;
    } else {
        lowJ = maxJ - 5;
    }
    if ((maxJ + 5) > w) {
        highJ = w;
    } else {
        highJ = maxJ + 5;
    } 

    writeln("lowI: ", lowI, " highI: ", highI, "lowJ: ", lowJ, " highJ: ", highJ);
    const subD: domain(2) = {lowI..highI, lowJ..highJ};
    for (i,j) in subD {
        sum -= result(i,j);
    }
    
    //Calculate average energy
    var energy : real;
    energy = sum/((h*w) - ((highI-lowI + 1)*(highJ-lowJ + 1)));
    // energy = sum/((h*w) - 121);
    writeln("energy: ", energy);
    var PCE = max/energy;
    writeln("pce: ", PCE);
    return PCE;
}
