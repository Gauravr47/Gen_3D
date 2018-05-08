#ifndef IMAGE_ARRAY
#define IMAGE_ARRAY
#include"CTimage.h"
#include"spline.h"


// array of our defined structure
  Image *imagearray; /* Imagearray is the structure that stores the new intermediate images */
  Image *inputarray; /* Input array is the struture which only stires the input images. The size of the array is equal to the no of input images*/
  int no_layers;
  int total_images;

  
// Structure to store traingles constructed usin marching cube


// functions on array of structure
 void preInitializeCheck(IPM*, char**);
 void allocateInputArray ( int , char**, int, int, double, double);
 void interpolate ( double, double, int);
 void cubicInterpolate(double, double, int)
 bool initialize_imagearray ( double, double, int, IPM&);
 void cropImagearray(int &, int &, IPM&)
 void deconstruct_arrays();
 void initialize_paddedarray(int, int , int, double ,double);
 void threshold(Image&, unsigned char);

#endif

