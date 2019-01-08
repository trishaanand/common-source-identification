void grayscale(int h, int w, double *image_grayscale, unsigned char *image_rgb) {
  for (int i = 0; i < h * w; i++) {
    double r = (double) image_rgb[i * 3 + 0];
    double g = (double) image_rgb[i * 3 + 1];
    double b = (double) image_rgb[i * 3 + 2];
    image_grayscale[i] = 0.299 * r + 0.587 * g + 0.114 * b;
  }
}
