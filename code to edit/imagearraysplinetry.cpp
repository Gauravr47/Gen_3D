

#include<sys/stat.h>
#include<iostream>
#include"imagearray.h"
#include"spline.h"
using namespace tk;
using namespace std;

/***************************************************************/
/*********** Allocate and deallocate arrays functions **********/
/***************************************************************/

/* Used to check the dimesnions of the input image. Initiates the image structure using check  */

void preInitializeCheck(int& columns, int& rows, int& maxVal,  char**argv)
{Image check;
 FILE *fcheck = fopen(argv[1],"rb");
 check.readHeader(fcheck);
 columns= check.numberOfColumns;
 rows= check.numberOfRows;
 maxVal= check.maxVal;
}

/*****************imagearray************************************/
/* Imagearray is the structure that stores the new intermediate images */

bool initialize_imagearray ( double y_scan, double  y_AM, int argc, int columns, int rows, int maxVal, double x_z_pixel)
{  /* Here y_AM gives us the  new y-resolution we expect from the 3-D printed parts, y_scan gives us the current y-resolution 
from the given ct images and total_images gives us the no of pgm images generated 
   */
   if (y_AM>y_scan)
  { cout<<"\n New resolution cannot be bigger than current resolution\n";
    return false; 
  }
   if (y_AM !=0)
   {  total_images = int((argc-2)*y_scan/y_AM)+1;
      imagearray =(Image*)malloc(sizeof(Image)*total_images);
      if(imagearray==NULL)
        { cout<<"\n Imagearray  allocation was not successfull \n";
        }
      for(int i=0; i<total_images;i++)
       { imagearray[i].initialize(columns, rows, maxVal, x_z_pixel, y_AM);
       }
   }
  else
   {  cout<<"\n New z-resolution cannot be zero";
      exit;
   }

  cout<<"\n Total no of images created are "<< total_images <<"\n";
  return true;
}

/************************ inputarray****************************/
/* Input array is the struture which only stires the input images. The size of the array is equal to the no of input images*/

void allocateInputArray (int argc, char**argv, int columns, int rows, int maxVal, double x_z_pixel,double y_pixel)
{ inputarray = (Image*) malloc (sizeof(Image)*(argc-1));

  for (int i=0; i <argc-1; i++)

        { inputarray[i].initialize(columns, rows, maxVal, x_z_pixel, y_pixel);
          inputarray[i].readImage (i+1, argv);
        }
}

/********************** paddedarray*****************************/
/* This array has a size of imagearray +2 in whhich the first and the last image has all the pixel defined as null pixel. Marching cube is only done on paddedarray. Here lenght denotes the original lenght of the array to be transfered to paddedarray */

int initialize_paddedarray(Image array[], int length)
{paddedarray =(Image*)malloc(sizeof(Image)*(length+2));
 int columns = array[0].numberOfColumns;
 int rows = array[0].numberOfRows;
 int maxVal = array[0].maxVal;
 double x_z_pixel = array[0].x_pixelSize;
 double y_pixel = array[0].y_pixelSize;
 for (int i =1; i< length+1; i++)
    {paddedarray[i].initialize(columns, rows, maxVal, x_z_pixel, y_pixel);
    for(int j=0; j<paddedarray[i].stringsize; j++)
{   paddedarray[i].imagedata[j]= array[i-1].imagedata[j];}
     paddedarray[i].paddImageBorders();
    }
 paddedarray[0].initialize(columns, rows, maxVal, x_z_pixel, y_pixel);
 paddedarray[length+1].initialize(columns, rows, maxVal, x_z_pixel, y_pixel);
 for(long i =0; i< paddedarray[0].stringsize; i++)
  {
   paddedarray[0].imagedata[i]=0;
   paddedarray[length+1].imagedata[i] = 0;
  }
 cout<<"\n Done initializing the paddedarray \n";
 return(length+2);
}

/* function to reallocate the paddedarray when new array is to shifted to paddedarray for mesh generation purposes */
 int  reallocate_paddedarray(Image array[], int length1, int length2)
{if (length1>length2)
  { for (int i = length2+2; i<length1; i++)
      { delete[] paddedarray[i].imagedata;
      }
  }
 paddedarray=(Image*)realloc(paddedarray, sizeof(Image)*(length2+2));
 int columns = array[0].numberOfColumns;
 int rows = array[0].numberOfRows;
 int maxVal = array[0].maxVal;
 double x_z_pixel = array[0].x_pixelSize;
 double y_pixel = array[0].y_pixelSize;
 for (int i =0; i<(length2+2); i++)
    { paddedarray[i].numberOfColumns =columns;
      paddedarray[i].numberOfRows =rows;
      paddedarray[i].maxVal= maxVal;
      paddedarray[i].x_pixelSize =x_z_pixel;
      paddedarray[i].z_pixelSize =x_z_pixel;
      paddedarray[i].y_pixelSize =y_pixel;
      paddedarray[i].reallocate(columns, rows);
      if(i!=0 && i!=(length2+1))
      { paddedarray[i].imagedata = array[i-1].imagedata;
        paddedarray[i].paddImageBorders(); 
      }
      if (i==0|| i==(length2+1))
      { for(long j =0; j< paddedarray[0].stringsize; j++)
           {
            paddedarray[i].imagedata[j]=0;
           }
      }
    }
  cout<<"\n exiting reallocation of paddedarray\n";
  return (length2+2);
}

/*********************Free all arrays **************************/

void  deconstruct_arrays (int argc, int length)
{  /*  free struct */
    for (int i=0; i<total_images; i++)
        { delete[] imagearray[i].imagedata;
        }
    for(int i=0; i<argc-1; i++)
        { delete[] inputarray[i].imagedata;
        }
/*    for(int i=0;i< length; i++)
        { delete[] paddedarray[i].imagedata;
        }*/
    free(imagearray);
    imagearray = NULL;
    cout<<"\n Imagearray  deallocation successfull \n";
    free( inputarray);
    inputarray = NULL;
    cout<<"\n Inputarray  deallocation successfull \n";
    free( paddedarray);
    paddedarray = NULL;
    cout<<"\n Paddedtarray  deallocation successfull \n";

}


/****************************************************************/
/************* Image arrays processing functions****************/
/***************************************************************/

/* Interpolate used for pixelwise linear interpolation of images*/

void interpolate ( double y_scan, double y_AM, int argc )
{ cout<<"\n Entered interpolation function \n";
  for( long i=0; i<(imagearray[0].stringsize);i++)
  {  unsigned char  pixeldifference =0;//difference in pixel value
     unsigned char basepixel =0;
     bool pixelpositive = true;
     int inputno = 0; // tracking position in inputarray
     double yabsolute = 0.0; // tracking spatial position of the image
     double yreference = y_scan; // reference used to keep a check on interpolation
     double ycorrection =0.0; // correction to yreference after it moves to next consecutive pair of images
     int inpcheck = 1.0;  // counter used to avoid calculating slope in every iteration

    for(int j=0; j<total_images; j++)
       { if(inputno <argc-1 && inpcheck!=inputno)
        {
          if(inputarray[inputno+1].imagedata[i]>inputarray[inputno].imagedata[i])
          { pixeldifference = inputarray[inputno+1].imagedata[i]-inputarray[inputno].imagedata[i];
            basepixel = inputarray[inputno].imagedata[i];
            pixelpositive= true;
          }
          if(inputarray[inputno+1].imagedata[i]<inputarray[inputno].imagedata[i])
          { pixeldifference = -inputarray[inputno+1].imagedata[i]+inputarray[inputno].imagedata[i];
            basepixel = inputarray[inputno].imagedata[i];
            pixelpositive= false;
          }
          if(inputarray[inputno+1].imagedata[i]==inputarray[inputno].imagedata[i])
          { pixeldifference = inputarray[inputno+1].imagedata[i]-inputarray[inputno].imagedata[i];
            basepixel = inputarray[inputno].imagedata[i];
            pixelpositive= true;
          }
	   inpcheck = inputno;
	}

	 if (yabsolute<= yreference && yreference <= (total_images+1)*y_AM)
	     {
              if(pixelpositive)
                { imagearray[j].imagedata[i]= basepixel +(yabsolute-ycorrection)*(pixeldifference)/y_scan;
		  yabsolute = yabsolute + y_AM;
                }

              if(!pixelpositive)
                { imagearray[j].imagedata[i]= basepixel -(yabsolute-ycorrection)*(pixeldifference)/y_scan;
                  yabsolute = yabsolute + y_AM;
                }
              }
	   if (yabsolute >yreference && yreference <= (total_images+1)*y_AM)
              { ycorrection = yreference;
		yreference = yreference +y_scan;
	        inputno++;
                
	      }
      }
   } cout<<"\n Exiting interpolate ";
}

/* Function to interpolate images using piiecwise cubic interpolation. Spline genrating function taken from tino klug repository. Here xin and yin are used to represent the spline cordinates in thier local system where x is the base axis  and y is the value axis...........
    y ^
      |    ''                 ''          '
      |   '  ''''         '  '  '''    '''
      |..'       ''.....'' '.      '..'
      |______________________________________________> 
                                                     X
Here xin dedones input grid points and xout denotes output points. Similarily yin and yout denotes corresponding values
 */

void cubicInterpolate (double y_scan, double y_AM, int argc)
{ cout<<"\n Entered cubic interpolation function \n";
  bool cubic_spline = true;
  for(long l =0; l<imagearray[0].stringsize;l++)
  { std:: vector<double> X(argc-1), Y(argc-1);
     double xin=0.0; // 
     double yin;
     double xout=0.0;
     double yout;
      for(int m =0; m<argc-1; m++)
     {
       X[m]= xin;
       xin = xin + y_scan;
       yin =(double) inputarray[m].imagedata[l];
       Y[m]= yin;
     }
    tk::spline s;
    s.set_points(X,Y, cubic_spline);
    for(int n =0; n<total_images; n++)
    {
      yout= s(xout);
      if (yout>255.0)
         { yout = 255.0;}
      if (yout < 0.0)
         { yout = 0.0;}
      imagearray[n].imagedata[l] = (unsigned char)yout;
      xout= xout+y_AM;
     }
   }
 cout <<"\n Done with cubic interpolation \n ";
}






/*********************marching cube ****************************/

void assignCube (GRIDCELL &grid,long iIndex, int yIndex )
{ double x_offset, y_offset, z_offset;
   x_offset= paddedarray[yIndex].x_pixelSize/2;
   z_offset = paddedarray[yIndex].z_pixelSize/2;
   y_offset = paddedarray[yIndex].y_pixelSize /2;

// Set virtual cube cordinates
   grid.p[0]= paddedarray[yIndex].stringToRealCords( iIndex, yIndex );
   grid.p[1]= paddedarray[yIndex].stringToRealCords( iIndex+1, yIndex );
   grid.p[2]= paddedarray[yIndex].stringToRealCords( iIndex+paddedarray[yIndex].numberOfColumns+1, yIndex );
     grid.p[3]= paddedarray[yIndex].stringToRealCords( iIndex+paddedarray[yIndex].numberOfColumns, yIndex );
   grid.p[4]= paddedarray[yIndex].stringToRealCords( iIndex, yIndex+1 );
    grid.p[5]= paddedarray[yIndex].stringToRealCords( iIndex+1, yIndex+1);
     grid.p[6]= paddedarray[yIndex].stringToRealCords( iIndex+paddedarray[yIndex].numberOfColumns+1, yIndex+1 );
    grid.p[7]= paddedarray[yIndex].stringToRealCords( iIndex+paddedarray[yIndex].numberOfColumns, yIndex+1 );
  
//Offset the virtual cube such that traingle base lies on the loaction of the slices
  for(int v =0; v<8; v++)
     { grid.p[v]= add (grid.p[v],x_offset, y_offset, z_offset);
     }

//Assign greyvalue to the cube vertices
  grid.val[0]= paddedarray[yIndex].imagedata[iIndex];
  grid.val[1]= paddedarray[yIndex].imagedata[iIndex+1];
  grid.val[2]= paddedarray[yIndex].imagedata[iIndex+paddedarray[yIndex].numberOfColumns+1];
  grid.val[3]= paddedarray[yIndex].imagedata[iIndex+paddedarray[yIndex].numberOfColumns];
  grid.val[4]= paddedarray[yIndex+1].imagedata[iIndex];
  grid.val[5]= paddedarray[yIndex+1].imagedata[iIndex+1];
  grid.val[6]= paddedarray[yIndex+1].imagedata[iIndex+paddedarray[yIndex].numberOfColumns+1];
  grid.val[7]= paddedarray[yIndex+1].imagedata[iIndex+paddedarray[yIndex].numberOfColumns];

}

/* Find the location of a traingle vertex on the edge of the virtual cube using interpolation */

Dp3 getEdgeIntersection(GRIDCELL &grid, int vertex1 , int vertex2, unsigned char Tlevel )
{  Dp3 v1,v2, e;
   v1 = grid.p[vertex1];
   v2= grid.p[vertex2];
   e.x= v1.x+(v2.x-v1.x)*(Tlevel - grid.val[vertex1])/(grid.val[vertex2] - grid.val[vertex1]);
   e.y=v1.y+(v2.y-v1.y)*(Tlevel - grid.val[vertex1])/(grid.val[vertex2] - grid.val[vertex1]);
   e.z= v1.z+(v2.z-v1.z)*(Tlevel - grid.val[vertex1])/(grid.val[vertex2] - grid.val[vertex1]);
   return(e);
}

//Get vertex from 
Dp3 getVertex(int edgeno, Dp3 * edgeIntersectList)
{ Dp3 e;
 e = edgeIntersectList[edgeno];
 return(e);
}
 
void  getNormal(TRIANGLE &tri)
{ Dp3 p1,p2,p3;
  p1 = tri.p[0];
  p2 = tri.p[1];
  p3 = tri.p[2];
  tri.normal.x = (p2.y -p1.y)*(p3.z -p1.z) - (p3.y -p1.y)*(p2.z -p1.z);
  tri.normal.y = (p2.z -p1.z)*(p3.x -p1.x) - (p2.x -p1.x)*(p3.z -p1.z);
  tri.normal.z = (p2.x -p1.x)*(p3.y -p1.y) - (p3.x -p1.x)*(p2.y -p1.y);
} 




/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/

long marchingCube(int length,unsigned char Tlevel)
{
   long nTriang=0;
   int cubeIndex, edgeIndex;
   Dp3 vertlist[12];
   Dp3 edgeIntersectList[12];
   GRIDCELL grid;


   extern int edgeTable[256];
   extern int triTable[256][16];
   extern int edgeConnection[12][12];


  for (int yIndex = 0; yIndex < length-1; yIndex++)
  {   int columns = paddedarray[yIndex].numberOfColumns;
      int rows = paddedarray[yIndex].numberOfRows;
      for ( long iIndex =0; iIndex< paddedarray[yIndex].stringsize-1; iIndex++)
       {  Ip2 check;
          check = paddedarray[yIndex].stringToCords( iIndex);
         if ((check.x < (columns-1)) && (check.z< (rows -1)))
          { assignCube(grid , iIndex, yIndex);
            //cout<< check.x << "\t" << check.z <<"\t"<<grid.p[0].x <<"\t" ;
    /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
    */
            cubeIndex = 0;
           for(int v = 0; v < 8; v++)
            {
                if(grid.val[v] > Tlevel) 
                     {  cubeIndex |= 1 <<v; }
            }
   /* Cube is entirely in/out of the surface */
           if (cubeIndex == 0)
             { continue;}
         //   cout << check.x <<"\t " << cubeIndex <<" \t" ;
            edgeIndex =edgeTable[cubeIndex]; 
   /* Find the vertices where the surface intersects the cube */
   //Then find the normal to the surface at those points
          for (int iEdge = 0; iEdge < 12; iEdge++)
            {
   //if there is an intersection on this edge
              if(edgeIndex & (1<<iEdge))
                {
                    edgeIntersectList[iEdge] = getEdgeIntersection(grid,edgeConnection[iEdge][0] , edgeConnection[iEdge][1], Tlevel );

                }
            }
//Generate triangles
           for(int iTriangle = 0; iTriangle < 5; iTriangle++)
           {
                TRIANGLE tri;
                if(triTable[cubeIndex][3*iTriangle] < 0)
                      {  break; }
                for(int iCorner = 0; iCorner < 3; iCorner++)
                   {
                        tri.p[iCorner] = getVertex(triTable[cubeIndex][3*iTriangle+iCorner],edgeIntersectList);
                   }
                getNormal( tri );
                triangles.push_back(tri);
                nTriang ++;
           }
        }
     }
  }

   return(nTriang);
 cout<<"\n Done Tessellating ";
}


/***************************************************************/
/**************** Image processing functions *******************/
/***************************************************************/

void smoothing(Image& img, double sigma,int kernalsize)
{ Image test;
  int columns=img.numberOfColumns;
  int rows = img.numberOfRows;
  int maxval = img.maxVal;
  int xgint;
  test.initialize(columns, rows, maxval);
  double r;
  double  s=2.0 * sigma * sigma;
  double sum =0.0;
  double xg =0.0;
  int width = test.numberOfColumns;
  double gKernal[kernalsize][kernalsize];
  double matrix[kernalsize][kernalsize];
  int k = kernalsize/2;
  Ip2 point,matrixmaker;
  long int  matrixindex;

 for(int x =-k; x<=k; x++)
  {  for(int z= -k; z<= k; z++)
     {  r= sqrt(x*x+z*z);
        gKernal[x+k][z+k]= (exp(-(r*r)/s))/(M_PI*s);
        sum+=gKernal[x+k][z+k];
     }
  }
 cout<<"\n Kernal sum is "<<sum<<"\n";
 for(int i =0; i<kernalsize; i++)
   { for(int j=0; j <kernalsize; j ++)
     {  gKernal[i][j]/=sum;
        cout<< gKernal[i][j]<<"\t";;
     } cout<<"\n";
   }


 for (long int i =0; i< test.stringsize-1; i++)
    {xg=0; 
     point = img.stringToCords(i);
     cout<<(int) img.imagedata[i];
     cout<< "\n";
     if(point.x <rows-k && point.z<columns-k && point.x>k && point.z>k)
       {   for(int x =-k; x<=k; x++)
             {  for(int z= -k; z<= k; z++)
                    { matrixmaker.x = point.x +x;
                      matrixmaker.z = point.z +z;
                      matrixindex =img.cordsToString(matrixmaker);
                      matrix[x][z] =img.imagedata[matrixindex];
                      cout<< matrix[x][z]<<"\t";
		      xg += gKernal[x][z]*matrix[x][z];
                    }
              }
       /*   for(int x =-k; x<=k; x++)
              {  for(int z= -k; z<= k; z++)
                    {  xg += gKernal[x][z]*matrix[x][z];
                    }
              } */

        }
        xgint = xg;
        test.imagedata[i] =(unsigned char) xgint;
        cout<<"\t\t"<<xg <<"\t"<<xgint <<"\n";
        xg =0;
    }

  for (long int i=0; i< test.stringsize; i++)
     { img.imagedata[i] = test.imagedata[i];
     }
 }



//Threshold on original image

void threshold( Image& img , unsigned char tLevel)
{ 
   for (long int j =0; j< img.stringsize; j++)
     {  if ((tLevel+(unsigned char)5) <= img.imagedata[j] || (tLevel-(unsigned char)5) <= img.imagedata[j] )
	  { img.imagedata[j] =(unsigned char) 255;
	  }
        else
	   { img.imagedata[j]=(unsigned char) 0;
	   }
     }
}

// Threshold and generate new image
Image thresholdseries (Image& img, unsigned char tLevel)
{ Image filt;
  filt.initialize(img.numberOfColumns, img.numberOfRows, img.maxVal);
 for (long int j =0; j< img.stringsize; j++)
     {  if (tLevel <= img.imagedata[j])
          { filt.imagedata[j] =(unsigned char) 255;
          }
        else
           { filt.imagedata[j]=(unsigned char) 0;
           }
     }
  return(filt);
}

// Cropping the image

void imageCrop (Image& img,int rows, int columns,int x_offset, int z_offset)
{ Image temp;
  temp.initialize(columns, rows, 255);
  Ip2 check, doublecheck;
  Ip2 startpoint, endpoint;
  char outimgname[20];
  long counter =0;
  long startindex;
  long endindex;
  startpoint.x= x_offset;
  startpoint.z= z_offset;
  startindex = img.cordsToString(startpoint);
  endpoint.x= startpoint.x+columns-1;
  endpoint.z= startpoint.z+rows-1;
  endindex = img.cordsToString(endpoint); 

  for(long i=startindex; i<=endindex; i++)
    {
       check = img.stringToCords(i);
     if( check.x >=startpoint.x && check.x <= endpoint.x && check.z >= startpoint.z  && check.z <= endpoint.z)
       {
         temp.imagedata[counter] = img.imagedata[i];
         counter++;
         doublecheck.x = check.x; doublecheck.z=check.z;
       }
   }
  cout<< counter;
   img.reallocate(columns, rows);
  for( long i = 0; i< img.stringsize; i++)
      {  img.imagedata[i]= temp.imagedata[i];
      }

 cout<< "\n Done cropping the images ";
}

// Invert the image
void inversion(Image& temp)
{ unsigned char val; 
  for (long int j =0; j<temp.stringsize; j++)
     {  val = temp.imagedata[j];
        temp.imagedata[j]= temp.maxVal -val;
     }
}

/***************************************************************/
/********************* write stl functions *********************/
/***************************************************************/

//ASCII STL
void writeSTL(char outname[])
{ FILE *fstl;
  fstl = fopen(outname, "w");
  fprintf( fstl, " solid ctscan \n");
  for (std::list<TRIANGLE> ::iterator it = triangles.begin(); it != triangles.end(); ++it)
   { TRIANGLE ves = *it;
     fprintf (fstl, " facet normal %e %e %e \n", ves.normal.x, ves.normal.y, ves.normal.z );
     fprintf ( fstl,"   outer loop \n" );
     fprintf ( fstl,"      vertex %e %e %e\n", ves.p[0].x, ves.p[0].y, ves.p[0].z );
     fprintf ( fstl,"      vertex %e %e %e\n", ves.p[1].x, ves.p[1].y, ves.p[1].z );
     fprintf ( fstl,"      vertex %e %e %e\n", ves.p[2].x, ves.p[2].y, ves.p[2].z ); 
     fprintf (fstl, "    endloop\n");
     fprintf (fstl, " endfacet\n");
  }
  fprintf(fstl,  " endsolid ctscan");
  fclose(fstl);
}

//Binary STL
void writeSTLbin(char outname[], long nTriang)
{ char head[80] =" GR ctscan";
  short attribute =0;
  float convertor;
  unsigned long n= (unsigned long)nTriang;
  FILE *fstl;
  fstl = fopen(outname, "wb");
  /*for(list<TRIANGLE>::iterator i = triangles.begin(); i != triangles.end(); ++i)
   {n=triangles.size();} */
  fwrite((char*) &head, sizeof(char),80,fstl);
  fwrite((char*) &n, sizeof(n),1,fstl);
  for (std::list<TRIANGLE> ::iterator it = triangles.begin(); it != triangles.end(); ++it)
   { TRIANGLE ves = *it;
     convertor = (float) ves.normal.x;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
    // fprintf(fstl, "%f", convertor);
     convertor = (float) ves.normal.y;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
  //   fprintf(fstl, "%f", convertor);
     convertor = (float) ves.normal.z;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[0].x;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
    // fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[0].y;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[0].z;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
    // fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[1].x;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[1].y;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[1].z;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[2].x;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[2].y;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
     convertor = (float) ves.p[2].z;
 fwrite (( char*) &convertor, sizeof(convertor),1, fstl);
   //  fprintf(fstl, "%f", convertor);
 fwrite (( char*) &attribute, sizeof(attribute),1, fstl);
 //    fprintf(fstl, "%h", attribute);
  }
  fclose(fstl);
}







/***************************************************************/
/*********************** main **********************************/
/***************************************************************/

 int main( int argc, char** argv)
{  bool structarray;
   unsigned char  Tlevel=0;
   int t,paddedarraylength, kernalsize;
   short crop,smooth, invert;
   double y_scan, y_AM, x_z_pixel, sigma;
   char outimgname[256];
   long nTriang;
   int columns, rows, maxVal,xCropOffset, zCropOffset;
   cout <<"\n Enter the current pixel size in mm \n";
   cin>>x_z_pixel;
   cout<<"\n Enter the current y-resolution in mm \n";
   cin>>y_scan;
   cout<<"\n Enter the desired y-resolution in mm\n";
   cin>>y_AM;
   cout<<"\n Enter Threshold value \n";
   cin >> t;
   Tlevel = (unsigned char)t;
   cout <<(int) Tlevel;
   preInitializeCheck(columns, rows, maxVal, argv);

   allocateInputArray ( argc, argv, columns, rows, maxVal, x_z_pixel,y_scan);


/*   cout<<"\n Do wou want to smooth the images ? \n Enter 1 for yes \n 2 for No \n";
   cin>>smooth;
   if (smooth ==1)
   {  cout <<"\n Enter kernal size (Must be an odd number) \n";
      cin>> kernalsize;
      cout<<"\n Enter sigma value \n";
      cin>> sigma;
      for(int i=0; i<argc-1; i++)
          {smoothing(inputarray[i],sigma, kernalsize);
          }
   } */

   for(int i=0; i<argc-1;i++)
      { snprintf(outimgname, sizeof(outimgname), "output/input_%d.pgm", i+1001);
        inputarray[i].writeImage(outimgname);
      }
    
   cout<<" \n Do you want to crop image? \n Enter \n 1 for YES \n 2 for NO \n";
   cin>>crop;
   if( crop == 1)
   { cout<< "\n Enter the desired number of rows \n" ;
     cin>> rows;
     cout<<"\n Enter the desired number of columns \n";
     cin >> columns;
     cout<<"\n Enter the offset in x and z direction respectively \n";
    cin >> xCropOffset >> zCropOffset;
    for(int i=0; i<argc-1; i++)
    {imageCrop(inputarray[i],rows, columns, xCropOffset, zCropOffset);
    }
   }

   cout<<"\n Do you want to invert the colors of the image \n  Enter \n 1 for YES \n 2 for NO \n";
   cin>> invert;
   if (invert==1)
    { for(int i=0; i<argc-1;i++)
      { inversion(inputarray[i]);
      }
    }
 
   structarray = initialize_imagearray( y_scan, y_AM, argc , columns, rows, maxVal, x_z_pixel);
   if(!structarray){ cout<<" \n Structure array was not initailized properly \n "; }
  else 
  { /*Switch between linear and Cubic by commenting one and uncommentig the other */
     interpolate( y_scan, y_AM, argc);
  // cubicInterpolate( y_scan, y_AM, argc); 

  for(int i=0; i<total_images;i++)
      { snprintf(outimgname, sizeof(outimgname), "output/Output_%d.pgm", i+1001);
        imagearray[i].writeImage(outimgname);
      }

/*  //threshold switch 
    for(int i=0; i<total_images;i++)
      { threshold(imagearray[i], Tlevel);
      }
*/

 // Stl of imagearray (output)
  paddedarraylength =initialize_paddedarray(imagearray,total_images);
  for(int i=0; i<paddedarraylength;i++)
      {threshold(paddedarray[i] , Tlevel);
        snprintf(outimgname, sizeof(outimgname), "output/padded1_%d.pgm", i+1001);
        paddedarray[i].writeImage(outimgname);
     }

  nTriang = marchingCube (paddedarraylength,Tlevel);
  snprintf(outimgname, sizeof(outimgname), "output/outstlbin.stl");
  writeSTLbin(outimgname, nTriang);
  cout <<" #triangles in output = "<< nTriang <<"\n ";

  deconstruct_arrays(argc, paddedarraylength);
  std::cout<<"\n Operation completed successfull ";
  return 0;
}
}



/****************************************************************/ 
/**************  tables to be used in functions*****************/ 
/***************************************************************/

 int edgeConnection[12][12] = 
{
        {0,1}, {1,2}, {2,3}, {3,0},
        {4,5}, {5,6}, {6,7}, {7,4},
        {0,4}, {1,5}, {2,6}, {3,7}
};
 int edgeTable[256]={
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

 int triTable[256][16] =
{       {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
        {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
        {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
        {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
        {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
        {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
        {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
        {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
        {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
        {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
        {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
        {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
        {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
        {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
        {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
        {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
        {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
        {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
        {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
        {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
        {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
        {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
        {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
        {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
        {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
        {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
        {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
        {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
        {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
        {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
        {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
        {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
        {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
        {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
        {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
        {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
        {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
        {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
        {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
        {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
        {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
        {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
        {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
        {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
        {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
        {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
        {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
        {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
        {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
        {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
        {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
        {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
        {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
        {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
        {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
        {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
        {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
        {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
        {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
        {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
        {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
        {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
        {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
        {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
        {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
        {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
        {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
        {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
        {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
        {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
        {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
        {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
        {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
        {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
        {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
        {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
        {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
        {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
        {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
        {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
        {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
        {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
        {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
        {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
        {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
        {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
        {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
        {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
        {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
        {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
        {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
        {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
        {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
        {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
        {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
        {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
        {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};




