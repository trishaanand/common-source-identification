/* Use the system C types to conveniently convert the C types to Chapel types 
 * and vice versa.
 */
use SysCTypes;

/* Use the read_jpg library for the RGB record.*/
use read_jpg;

/* Instruct the compiler to compile the following files and link them in the 
 * binary.  This C-code will apply filters to compute a PRNU noise pattern.
 */
require "../c/prnu.c", "../c/prnu.h", 
  "../c/fastnoise.c", "../c/fastnoise.h",
  "../c/grayscale.c", "../c/grayscale.h",
  "../c/util.c", "../c/util.h",
  "../c/wiener.c", "../c/wiener.h",
  "../c/zeromean.c", "../c/zeromean.h";

/* Declare a Chapel record to mirror a C struct for holding temporary data.
 *
 * To apply all filters to obtain a PRNU noise pattern, we need temporary data.
 * This record keeps track of pointers to the temporary data and the FFT plans.
 * 
 * We can initialize and destroy this temporary data with the pointers using the
 * initialize and destroy functions below.
 */
extern record prnu_data {
  var h : c_int;
  var w : c_int;
  var h_w_double1 : c_ptr(c_double);
  var h_w_double2 : c_ptr(c_double);
  var h_w_double3 : c_ptr(c_double);
  var h_w_double_border : c_ptr(c_double);
  // The following two declarations are actually: *fftwf_fcomplex
  var h_w_forward : c_ptr(c_double); 
  var h_w_backward : c_ptr(c_double);
  // The following two declaration are actually: fftwf_plan, which are 
  // pointers to a struct
  var plan_forward : c_void_ptr;
  var plan_backward : c_void_ptr;
}

/* Declare three functions to interface with the C code.
 * We will create more Chapel-like wrappers below.
 */
extern proc prnu_initialize(h : c_int, w : c_int, 
			    ref prnu_data_ref : prnu_data) : void;

extern proc prnu_destroy(ref prnu_data_ref : prnu_data) : void;

extern proc prnu_execute(ref prnu : c_double, ref image_rgb : c_uchar, 
			 ref prnu_data_ref : prnu_data) : void;

/* Initialize the prnu_data data structure */
proc prnuInit(h : int, w : int, ref data : prnu_data) {
  prnu_initialize(h : c_int, w : c_int, data);
}

/* Compute a PRNU noise pattern from an image with RGB data. 
 *
 * This function assumes that prnu and image are arrays with domains that are
 * sufficient to store all the data.  Concretely this means that the domain
 * should have h*w reals and RGB pixels respectively and the C function
 * assumes that it is contiguous memory.  Note that h and w are defined in the
 * prnu_data record.
 */
proc prnuExecute(prnu : [] real, image : [] RGB, ref data : prnu_data) {
  prnu_execute(prnu[prnu.domain.low], image[image.domain.low].r, data);
}

/* Destroy the memory that is reserved for doing PRNU computations.
 *
 * The memory that has been allocated with prnuInit() will be destroyed using
 * this function.
 */
proc prnuDestroy(ref data : prnu_data) {
  prnu_destroy(data);
}

