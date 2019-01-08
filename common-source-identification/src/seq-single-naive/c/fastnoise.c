#include <math.h>

#include "util.h"

#define EPS 1.0

/*
 * Vertically computes a local gradient for each pixel in an image.  
 * 
 * Takes forward differences for first and last row. Takes centered differences
 * for interior points.
 */
void convolveVertically(int h, int w, double *output, double *input) {
  for (int j = 0; j < w; j++) {
    output[0 * w + j] += input[1 * w + j] - input[0 * w + j];
    output[(h - 1) * w + j] += input[(h - 1) * w + j] - input[(h - 2) * w + j];
    
    for (int i = 1; i < h - 1; i++) {
      output[i * w + j] += 0.5f * (input[(i + 1) * w + j] - input[(i - 1) * w + j]);
    }
  }
}

/*
 * Horizontally computes a local gradient for each pixel in an image. 
 * 
 * Takes forward differences for first and last element.  Takes centered
 * differences for interior points.
 */
void convolveHorizontally(int h, int w, double *output, double *input) {
  for (int i = 0; i < h; i++) {
    output[i * w + 0] += input[i * w + 1] - input[i * w + 0];
    output[i * w + w - 1] += input[i * w + w - 1] - input[i * w + w - 2];

    for (int j = 1; j < w - 1; j++) {
      output[i * w + j] += 0.5f * (input[i * w + j + 1] - input[i * w + j - 1]);
    }
  }
}

/*
 * Normalizes gradient values in place.
 */
void normalize(int h, int w, double *dxs, double *dys) {
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      double dx = dxs[i * w + j];
      double dy = dys[i * w + j];

      double norm = (double) sqrtf((dx * dx) + (dy * dy));
      double scale = 1.0f / (EPS + norm);

      dxs[i * w + j] = scale * dx;
      dys[i * w + j] = scale * dy;
    }
  }
}

void fastnoise(int h, int w, double *input_output, double *dxs, double *dys) {
  zero(h * w, dxs);
  zero(h * w, dys);

  convolveHorizontally(h, w, dxs, input_output);
  convolveVertically(h, w, dys, input_output);
  normalize(h, w, dxs, dys);

  zero(h * w, input_output);
  
  convolveHorizontally(h, w, input_output, dxs);
  convolveVertically(h, w, input_output, dys);
}
