//Header file for Image structure

#ifndef _CT_IMAGE
#define _CT_IMAGE

#include"points.h"
#include<iostream>
#include<cstdlib>
#include<cstdio>
//#include<math.h>
//#include <fstream> 
using namespace std;

struct ImagePrivateMembers {
	int columns = 0;
	int rows = 0;
	int maxVal = 0;
	double x_z_pixelSize = 0.0;
	double y_pixelSize = 0.0;

} typedef IPM;

// Structure for storing the image
 class Image
{ 
   int numberOfColumns, numberOfRows, maxVal;
//variables for store spatial resolution in all directions
   double  x_z_pixelSize, y_pixelSize; 
// String to store image header for Image writting functions
   char * imageheader;
//Variable used to track the position of the image
   int z_level;
// Variable storing imagedata size (used in intializing the array
   long stringsize;
   //functions
// Constructers 3 different types
   

public:
  // Array to store pixel data
  unsigned char * imagedata;
  Image() {
	  numberOfColumns = 0;
	  numberOfRows =0;
	  x_z_pixelSize =1;
	  maxVal =255;
          };

  void initialize(const IPM &);
  //Accesors
  int getColumns() { return numberOfColumns; };
  int getRows() { return numberOfRows; };
  int getMaxVal() { return maxVal; };
  long getStringSize() { return stringSize; };

//Function to reallocate the image for cropping or maginification purposes
   void reallocate(const int&,const int&);
//Image reading and writing functions
   void readImage (int ,char**);
   void readHeader(FILE*);
   void writeImage (char *);
   void writeHeader();
//Function to calculate x and y pixel position (in plane) from array index number
   Ip2  stringToCords (const long int&);
//Function to calculate array index number from x and y pixel position
   long int cordsToString (Ip2&);
//Function to calculate pixel cordinates in global 3D cartesian system
   Dp3  stringToRealCords (const long int&, const int&);
//Function to calculate array index from 3D cartesian cordinates
   int realCordsToString (double, double, char*);
//Function to padd croped images with a null value border
   void paddImageBorders();
//Destructer
   ~Image(){};
};


#endif
