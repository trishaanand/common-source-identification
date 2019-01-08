#include <stdlib.h>

#include <fftw3.h>

#include "prnu.h"
#include "grayscale.h"
#include "fastnoise.h"
#include "zeromean.h"
#include "wiener.h"

double *alloc_2d_double_array(int h, int w) {
  return (double *) malloc(h * w * sizeof(double));
}

void prnu_initialize(int h, int w, prnu_data *prnu_data) {
  prnu_data->h = h;
  prnu_data->w = w;
  prnu_data->h_w_double1 = alloc_2d_double_array(h, w);
  prnu_data->h_w_double2 = alloc_2d_double_array(h, w);
  prnu_data->h_w_double3 = alloc_2d_double_array(h, w);
  int border_size = MAX_FILTER_SIZE / 2;
  prnu_data->h_w_double_border = alloc_2d_double_array(h + 2 * border_size, 
						     w + 2 * border_size);
  fftw_complex *in_out = 
    (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * h * w);
  prnu_data->h_w_forward = in_out;
  prnu_data->plan_forward =
    fftw_plan_dft_2d(h, w, in_out, in_out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  in_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * h * w);
  prnu_data->h_w_backward = in_out;
  prnu_data->plan_backward =
    fftw_plan_dft_2d(h, w, in_out, in_out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void prnu_destroy(prnu_data *prnu_data) {
  fftw_destroy_plan(prnu_data->plan_backward);
  fftw_destroy_plan(prnu_data->plan_forward);
  fftw_free(prnu_data->h_w_backward);
  fftw_free(prnu_data->h_w_forward);
  free(prnu_data->h_w_double_border);
  free(prnu_data->h_w_double3);
  free(prnu_data->h_w_double2);
  free(prnu_data->h_w_double1);
}

void prnu_execute(double *prnu, unsigned char *image_rgb, prnu_data *prnu_data) {
  int h = prnu_data->h;
  int w = prnu_data->w;
  grayscale(h, w, prnu, image_rgb);
  fastnoise(h, w, prnu, prnu_data->h_w_double1, prnu_data->h_w_double2);
  zeromean(h, w, prnu, prnu_data->h_w_double1);
  wiener(h, w, prnu, prnu_data->h_w_forward,
	 prnu_data->h_w_double1, prnu_data->h_w_double2, 
	 prnu_data->h_w_double_border, prnu_data->h_w_double3,
	 prnu_data->h_w_backward,
	 prnu_data->plan_forward, prnu_data->plan_backward);
}
