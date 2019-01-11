use FFTW_MT;
use Math;
use Time;

/* prnu is passed by reference. Hence not returning any value from the FFT calculation */
proc calculateFFT(prnu : [] complex, sign : c_int) {
    var fftPlan = plan_dft(prnu, prnu, sign, FFTW_ESTIMATE);
    execute(fftPlan);
}

proc computeEverything(h : int, w : int, prnuComplex : [] complex, prnuRotComplex : [] complex) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    var resultComplex : [imageDomain] complex;
    var result : [imageDomain] real;

    var t1Timer, t2Timer, t3Timer, t4Timer, t5Timer  : Timer;

    t1Timer.start();
    // Calculate the point wise product of both matrices
    for (i,j) in imageDomain do {
        resultComplex(i,j) = prnuComplex(i,j) * prnuRotComplex(i,j);
    }
    t1Timer.stop();

    t2Timer.start();
    // Calculate inverse FFT of the result
    calculateFFT(resultComplex, FFTW_BACKWARD);
    t2Timer.stop();

    // Scale the result
    t3Timer.start();
    for(i,j) in imageDomain {
        result(i,j) = resultComplex(i,j).re;
        result(i,j) = ( result(i,j) * result(i,j) )/ ((h*w) * (h*w));
    }
    t3Timer.stop();
    
    
    t4Timer.start();
    var max, sum : real;
    sum = 0.0;
    var maxI, maxJ : int;
    for (i,j) in imageDomain do {
        //Find the peak
        if (max < result(i,j)) {
            max = result(i,j);
            maxI = i;
            maxJ = j;
        }
    }
    t4Timer.stop();
    
    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w) then w else maxJ + 5;

        
    t5Timer.start();
    for (i,j) in imageDomain do {
        if(!(i > lowI && i < highI && j > lowJ && j < highJ )) {
            sum += result(i,j);
        }
    }
    t5Timer.stop();

    //Calculate average energy
    var energy : real;
    // energy = sum/((h*w) - ((highI-lowI + 1)*(highJ-lowJ + 1)));
    energy = sum/((h*w) - 121);
    var PCE = (max * max) / energy;
    return (PCE, t1Timer.elapsed(), t2Timer.elapsed(), t3Timer.elapsed(), t4Timer.elapsed(), t5Timer.elapsed());
}
