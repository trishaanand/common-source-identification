#ifndef WIENER_H
#define WIENER_H

#define NR_FILTER_SIZES 4
#define MAX_FILTER_SIZE 9

#include <fftw3.h>

/*
 * Applies the Wiener Filter in-place. 
 * 
 * This function needs several temporary data.  Note that
 * temp_squared_magnitudes_border should be an array with a size of h and w
 * where the width of the maximum border applied to all 4 sides.  It can then
 * be reused for convolutions with smaller border sizes.
 */
void wiener(int h, int w, double *input_output,
	    fftw_complex *temp_forward,
	    double *temp_squared_magnitudes,
	    double *temp_variance_estimates, 
	    double *temp_squared_magnitudes_border,
	    double *temp_convolution,
	    fftw_complex *temp_backward,
	    fftw_plan plan_forward, fftw_plan plan_backward);

#endif
