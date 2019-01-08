#ifndef UTIL_H
#define UTIL_H

#include <fftw3.h>

/*
 * Initialize array to a value.
 */
void initialize_array(int n, double *array, double value);

/*
 * Zero an array.
 */
void zero(int n, double *array);

/*
 * Transpose input to output.
 *
 * The dimensions of input are assumed to be h,w, the dimensions of input will
 * be w,h.
 */
void transpose(int h, int w, double *output, double *input);

/*
 * Apply convolution.
 *
 * It assumes that output is of size h * w, while the input has dimensions: h +
 * 2 * filter_size/2, w + 2 * filter_size/2.
 */
void convolve(int h, int w, int filter_size, double *output, double *input);

/*
 * Compute the pointwise minimum of input_output and input.
 *
 * The result is stored in input_output.
 */
void minimum_array(int n, double *input_output, double *input);

/*
 * Make complex numbers from reals.
 *
 * This function adds a zero for the imaginary part.
 */
void to_complex(int n, fftw_complex *output, double *input);

/*
 * Make real numbers from complex ones.
 *
 * This function takes the real part of the complex numbers.
 */
void to_real(int n, double *output, fftw_complex *input);

/*
 * Sum an array to a single value.
 */
double sum_array(int n, double *input);

/*
 * Perform a pointwise multliply of two arrays.
 *
 * The result is stored in input_output.
 */
void multiply_array(int n, double *input_output, double *input);

/*
 * Copy input to output, taking into account a border.
 *
 * The rest of the input will be initialized to 0.0.
 */
void copy_with_border(int h, int w, int size_border, double *output, double *input);

void write_2d_array(int h, int w, double *array, char *filename);

#endif
