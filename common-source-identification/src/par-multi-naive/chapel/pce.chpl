/* This function returns the max value and it's indices from a given 2D real matrix */
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

/* This function computes the PCE value from a given 2D complex matrix */
proc computePCE(h : int, w : int, resultComplex : [] complex) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    
    var result : [imageDomain] real = (resultComplex.re * resultComplex.re) / ((h*w) * (h*w));
    
    var (max, lowI, highI, lowJ, highJ) = getMax(result, h, w);
    var sum = + reduce result;
    
    var innerDomain : subdomain(imageDomain) = {lowI..highI-1, lowJ..highJ-1};
    var ignoreSum = + reduce result[innerDomain];
    
    sum = sum - ignoreSum;

    var PCE = max * ((h*w) - 121) / sum;
    return PCE;
}
