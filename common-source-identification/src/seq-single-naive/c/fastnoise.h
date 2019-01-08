#define EPS 1.0f

/* 
 * Apply the fastnoise filter.
 *
 * Note that this is an in-place filter that needs the two temporary buffers
 * dxs and dys that will be overwritten.
 */
void fastnoise(int h, int w, double *input_output, double *dxs, double *dys);
