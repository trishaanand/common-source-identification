/* Use the system C types to conveniently convert the C types to Chapel types 
 * and vice versa.
 */
use SysCTypes;

/* Instruct the compiler to compile the following files and link them in the 
 * binary.  This C-code will load JPG images.
 */
require "../c/read_jpg.c", "../c/read_jpg.h";

/* Declare two functions to interface with the C code.
 * We will create more Chapel-like wrappers below.
 */
extern proc get_dimensions_jpg(ref height: c_uint, ref width: c_uint, 
			       filename: c_string) : c_int;

extern proc read_jpg(ref buffer: c_uchar, filename : c_string) : c_int;


/* Represent a pixel in an image. */
record RGB {
  var r : c_uchar;
  var g : c_uchar;
  var b : c_uchar;
}

/* Retrieve the dimensions of a JPG file.
 * 
 * Given a JPG filename, return the dimensions of the file in a pair.
 * It will throw an error when a failure occurs.
 */
proc getDimensionsJPG(fileName : string) : 2*int throws {
  var h, w : c_uint;
  var result = get_dimensions_jpg(h, w, fileName.c_str());
  if (result != 0) {
    stderr.writeln("Cannot read JPG image ", fileName);
    throw new Error();
  }
  return (h:int, w:int);
}

/* Read in a JPG file into an image.
 * 
 * Given a JPG filename, read the image into image.
 * It will throw an error when a failure occurs.
 *
 * This function assumes that image is an array with a domain that is 
 * sufficient to store all the pixels (in RGB).  Concretely this means that the
 * domain should have h*w RGB pixels and the C function assumes that it is 
 * contiguous memory.
 */
proc readJPG(ref image : [] RGB, fileName : string) throws {
  ref startImage : c_uchar = image[image.domain.low].r;
  if (read_jpg(startImage, fileName.c_str()) != 0) {
    stderr.writeln("Cannot read JPG image ", fileName);
    throw new Error();
  }
}