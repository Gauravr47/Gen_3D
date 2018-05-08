#ifndef MARCH_CUBE
#define MARCH_CUBE

#include"imagearray.h"
#include <list>

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

long  marchingCube(Image[], int, unsigned char);
void assignCube(Image[], GRIDCELL &, long, int);
Dp3 getEdgeIntersection(GRIDCELL &, int, int);
Dp3 getVertex(int, Dp3*);
void getNormal(TRIANGLE &);
void writeSTL();
#endif // !MARCH_CUBE

