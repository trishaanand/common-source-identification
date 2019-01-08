/*
 * Applies the Zero Mean Total filter on the CPU. 
 * 
 * The filter is done in-place, but it nees a temporary buffer for the
 * transpose.
 */
void zeromean(int h, int w, double *input_output, double *transposed);
