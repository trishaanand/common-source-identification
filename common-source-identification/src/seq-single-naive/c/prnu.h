#ifndef PRNU_H
#define PRNU_H

#include <fftw3.h>

typedef struct prnu_data {
  int h;
  int w;
  double *h_w_double1;
  double *h_w_double2;
  double *h_w_double3;
  double *h_w_double_border;
  fftw_complex *h_w_forward;
  fftw_complex *h_w_backward;
  fftw_plan plan_forward;
  fftw_plan plan_backward;
} prnu_data;

void prnu_initialize(int h, int w, prnu_data *prnu_data);

void prnu_destroy(prnu_data *prnu_data);

void prnu_execute(double *prnu, unsigned char *image_rgb, prnu_data *prnu_data);

#endif
