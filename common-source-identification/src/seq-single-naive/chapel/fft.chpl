use FFTW_MT;

proc calculateFFT(h : int, w : int, prnu : [] real, prnuRot : [] real) {
    const imageDomain: domain(2) = {0..#h,0..#w};
    var prnuComplex, prnuRotComplex, resultComplex : [imageDomain] complex;

    for (i,j) in imageDomain do {
        prnuComplex(i,j) = prnu(i,j) + 0i;
        prnuRotComplex(i,j) = prnuRot(i,j) + 0i;
    }

    // TODO: These plans are expensive. Save it for them to be used later
    var forwardPrnu = plan_dft(prnuComplex, prnuComplex, FFTW_FORWARD, FFTW_ESTIMATE);
    var forwardPrnuRot = plan_dft(prnuRotComplex, prnuRotComplex, FFTW_FORWARD, FFTW_ESTIMATE);
    
    execute(forwardPrnu);
    execute(forwardPrnuRot);
    
    for (i,j) in imageDomain do {
        // resultComplex(i,j).re = 1.0;
        // resultComplex(i,j).im = 1.0;
    }
}