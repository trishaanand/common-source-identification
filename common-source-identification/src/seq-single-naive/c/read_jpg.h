#ifndef READ_JPG_H
#define READ_JPG_H

int get_dimensions_jpg(unsigned int *height, unsigned int *width, const char *filename);

int read_jpg (unsigned char *buffer, const char *filename);

#endif
