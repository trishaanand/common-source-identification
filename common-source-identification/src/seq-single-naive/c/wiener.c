#include <float.h>

#include <fftw3.h>

#include "util.h"

#include "wiener.h"

const int FILTER_SIZES[] = { 3, 5, 7, 9 };

/*
 * Computes the square of each frequency and stores the result as a real.
 */
void squared_magnitudes(int h, int w, double *output, 
			fftw_complex *frequencies) {
  for (int i = 0; i < h * w; i++) {
      double re = frequencies[i][0];
      double im = frequencies[i][1];
      output[i] = (re * re) + (im * im);
  }
}

/*
 * This function scales the frequencies in input with a combination of the global variance and an estimate for the local
 * variance at that position. Effectively this cleans the input pattern from low frequency noise.
 */
void scale_with_variances(int h, int w, fftw_complex *scaled, 
			  fftw_complex *input, double *variance_estimates, 
			  double variance) {
  
  for (int i = 0; i < h * w; i++) {
    double estimate = variance_estimates[i];
    double scale = variance / (variance > estimate ? variance : estimate);
    scaled[i][0] = input[i][0] * scale;
    scaled[i][1] = input[i][1] * scale;
  }
}

/*
 * Estimates the minimum local variances by applying all filters.
 *
 */
void variance_estimates(int h, int w, double *variance_estimates, 
			double *squared_magnitudes, 
			double *squared_magnitudes_with_border,
			double *temp_convolution) {
  
  initialize_array(h * w, variance_estimates, DBL_MAX);
  for (int i = 0; i < NR_FILTER_SIZES; i++) {
    int filter_size = FILTER_SIZES[i];
    copy_with_border(h, w, filter_size / 2, 
		     squared_magnitudes_with_border, squared_magnitudes);
    convolve(h, w, filter_size, temp_convolution, 
	     squared_magnitudes_with_border);
    minimum_array(h * w, variance_estimates, temp_convolution);
  }
}

void wiener(int h, int w, double *input_output,
	    fftw_complex *temp_forward,
	    double *temp_squared_magnitudes,
	    double *temp_variance_estimates, 
	    double *temp_squared_magnitudes_border,
	    double *temp_convolution,
	    fftw_complex *temp_backward,
	    fftw_plan plan_forward, fftw_plan plan_backward) {
  to_complex(h * w, temp_forward, input_output);

  fftw_execute(plan_forward);

  squared_magnitudes(h, w, temp_squared_magnitudes, temp_forward);

  variance_estimates(h, w, temp_variance_estimates, temp_squared_magnitudes,
  		     temp_squared_magnitudes_border, temp_convolution);

  int n = w * h;
  multiply_array(n, input_output, input_output);
  double variance = (sum_array(n, input_output) * n) / (n - 1);
  scale_with_variances(h, w, temp_backward, temp_forward,
  		       temp_variance_estimates, variance);

  fftw_execute(plan_backward);

  to_real(h * w, input_output, temp_backward);
}
