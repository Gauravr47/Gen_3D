#ifndef IMAGE_ARRAY
#define IMAGE_ARRAY
#include <list>
#include"CTimage.h"

// array of our defined structure
  Image *imagearray;
  Image *inputarray;
  Image *paddedarray;
  int no_layers;
  int total_images;


// Structure to store traingles constructed usin marching cube
typedef struct {
   Dp3 p[3];
   Dp3 normal;
} TRIANGLE;

//list of all traingles structures
std::list<TRIANGLE> triangles;

// Structure used as the virtual cube in the marching cube algorithm
typedef struct {
   Dp3 p[8];
   unsigned char val[8];
} GRIDCELL;


// Fuctions for mesh generation

long  marchingCube(Image[],int,unsigned char);
void assignCube (Image[],GRIDCELL &, long  , int);
Dp3 getEdgeIntersection(GRIDCELL &, int , int );
Dp3 getVertex(int, Dp3*);
void getNormal( TRIANGLE & );
void writeSTL();

// functions on array of structure
 void preInitializeCheck(int&, int&, int& , char**)
 void allocateInputArray ( int , char**, int, int, double, double);
 void interpolate ( double, double, int);
 bool initialize_imagearray ( double, double, int, int, int, int, double, double);
 void deconstruct_arrays();
 void initialize_paddedarray(int, int , int, double ,double);
 void threshold(Image&, unsigned char);

#endif

