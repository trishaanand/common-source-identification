#include <stdio.h>

#include "util.h"

void initialize_array(int n, double *array, double value) {
  for (int i = 0; i < n; i++) {
    array[i] = value;
  }
}

void zero(int n, double *array) {
  initialize_array(n, array, 0.0f);
}

void transpose(int h, int w, double *output, double *input) {
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      output[j * h + i] = input[i * w + j];
    }
  }
}


void convolve(int h, int w, int filter_size, double *output, double *input) {
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      double sum = 0.0f;
      for (int fi = 0; fi < filter_size; fi++) {
	for (int fj = 0; fj < filter_size; fj++) {
	  sum += input[(i + fi) * (w + (filter_size / 2) * 2) + (j + fj)];
	}
      }
      output[i * w + j] = sum / (filter_size * filter_size);
    }
  }
}

void minimum_array(int n, double *input_output, double *input) {
  for (int i = 0; i < n; i++) {
    double a = input_output[i];
    double b = input[i];
    input_output[i] = a < b ? a : b;
  }
}

void to_complex(int n, fftw_complex *output, double *input) {
  for (int i = 0; i < n; i++) {
    output[i][0] = input[i];
    output[i][1] = 0;
  }
}

void to_real(int n, double *output, fftw_complex *input) {
  for (int i = 0; i < n; i++) {
    output[i] = input[i][0] / n;
  }
}

double sum_array(int n, double *input) {
  double sum = 0.0f;

  for (int i = 0; i < n; i++) {
    sum += input[i];
  }

  return sum;
}

void multiply_array(int n, double *input_output, double *input) {
  for (int i = 0; i < n; i++) {
    input_output[i] *= input[i];
  }
}

void copy_with_border(int h, int w, int size_border, 
		      double *output, double *input) {
  // h and w with border size
  int h_b = h + 2 * size_border;
  int w_b = w + 2 * size_border;
  initialize_array(h_b * w_b, output, 0.0);
  
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      output[(i + size_border) * w_b + (j + size_border)] = input[i * w + j];
    }
  }
}

void write_2d_array(int h, int w, double *array, char *filename) {
  FILE *fp;
  fp=fopen(filename, "w+");

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      fprintf(fp, "%f ", array[i * w + j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}
