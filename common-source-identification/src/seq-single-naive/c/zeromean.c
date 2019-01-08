#include "util.h"

/*
 * Whether something is divisble by two.
 */
int is_divisible_by_two(int value) {
  return (value & 1) == 0;
}


/*
 * Applies an in-place zero mean filtering operation to each column. 
 * 
 * First two mean values are computed, one for even and one for odd
 * elements, for each column in the image. Then, the corresponding mean
 * value is subtracted from each pixel value in the image.
 */
void mean_vertically(int h, int w, double* input_output) {
  for (int j = 0; j < w; j++) {
    double sum_even = 0.0f;
    double sum_odd = 0.0f;

    for (int i = 0; i < h - 1; i += 2) {
      sum_even += input_output[i * w + j];
      sum_odd += input_output[(i + 1) * w + j];
    }
    if (!is_divisible_by_two(h)) {
      sum_even += input_output[(h - 1) * w + j];
    }

    double mean_even = sum_even / ((h + 1) / 2);
    double mean_odd = sum_odd / (h / 2);

    for (int i = 0; i < h - 1; i += 2) {
      input_output[i * w + j] -= mean_even;
      input_output[(i + 1) * w + j] -= mean_odd;
    }
    if (!is_divisible_by_two(h)) {
      input_output[(h - 1) * w + j] -= mean_even;
    }
  }
}

void zeromean(int h, int w, double *input_output, double *transposed) {
  mean_vertically(h, w, input_output);
  transpose(h, w, transposed, input_output);
  mean_vertically(w, h, transposed);
  transpose(w, h, input_output, transposed);
}
