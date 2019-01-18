use FFTW;
use Math;
use Time;

proc computePCE(h : int, w : int, resultComplex : [] complex) {
    const imageDomain: domain(2) = {0..#h, 0..#w};
    var result : [imageDomain] real;

    // Scale the result
    for(i,j) in imageDomain {
        result(i,j) = ( resultComplex(i,j).re * resultComplex(i,j).re )/ ((h*w) * (h*w));
    }
    
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
    
    //In the result matrix, remove 11x11 elements around the max from the total sum
    var lowI, highI, lowJ, highJ : int;
    lowI = if((maxI-5) < 0) then 0 else maxI -5 ;
    highI = if ((maxI + 5) > h) then h else maxI + 5;
    lowJ = if ((maxJ-5) < 0) then 0 else maxJ -5;
    highJ = if ((maxJ + 5) > w) then w else maxJ + 5;

    for (i,j) in imageDomain do {
        if(!(i > lowI && i < highI && j > lowJ && j < highJ )) {
            sum += result(i,j);
        }
    }

    //Calculate average energy
    var energy : real;
    energy = sum/((h*w) - 121);
    var PCE = max / energy;
    return PCE;
}
