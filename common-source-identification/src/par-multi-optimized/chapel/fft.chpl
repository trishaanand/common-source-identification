use FFTW;
use Time;

proc getMax(result : [] real, h : int, w : int) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    var lowI, highI, lowJ, highJ : int;
    var (max, maxLoc) = maxloc reduce zip(result, imageDomain);

    var maxI = maxLoc(1);
    var maxJ = maxLoc(2);
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w) then w else maxJ + 5;
    
    return (max, lowI, highI, lowJ, highJ);
}

proc computePCE(h : int, w : int, resultComplex : [] complex, ref pce : real) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    
    var result : [imageDomain] real;
    forall (i,j) in imageDomain {
        result(i,j) = (resultComplex(i,j).re * resultComplex(i,j).re) / ((h*w) * (h*w));
    }
    
    var (max, lowI, highI, lowJ, highJ) = getMax(result, h, w);
    var sum, ignoreSum : real;
    forall (i,j) in imageDomain with (+ reduce sum) {
        sum += result(i,j);
    }
    
    var innerDomain : subdomain(imageDomain) = {lowI..highI-1, lowJ..highJ-1};
    
    forall (i,j) in innerDomain with (+ reduce ignoreSum) {
        ignoreSum += result(i,j);
    }
    
    pce = max * ((h*w) - 121) / (sum - ignoreSum);
}