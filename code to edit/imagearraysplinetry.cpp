

#include<sys/stat.h>
#include"imagearray.h"
#include<cstdio>
using namespace std;

//using namespace tk;
/***************************************************************/
/*********** Allocate and deallocate arrays functions **********/
/***************************************************************/

/* Used to check the dimesnions of the input image. Initiates the image structure using check  */

void preInitializeCheck(IPM * ipm,  char**argv)
{Image check;
 FILE *fcheck = fopen(argv[1],"rb");
 check.readHeader(fcheck);
 ipm->columns= check.getColumns();
 ipm->rows= check.getRows();
 ipm->maxVal= check.getMaxVal();
}

/*****************imagearray************************************/
/* Imagearray is the structure that stores the new intermediate images */

bool initialize_imagearray (const double& y_scan, int & argc, IPM & ipm)
{  /* Here y_AM gives us the  new y-resolution we expect from the 3-D printed parts, y_scan gives us the current y-resolution 
from the given ct images and total_images gives us the no of pgm images generated 
   */
   if (ipm.y_pixelSize >y_scan)
  { puts("\n New resolution cannot be bigger than current resolution\n");
    return false; 
  }
   if (ipm.y_pixelSize !=0)
   {  total_images = int((argc-2)*y_scan/ipm.y_pixelSize)+1+2;
      imagearray =(Image*)malloc(sizeof(Image)*total_images);
      if(imagearray==NULL)
        { puts("\n Imagearray  allocation was not successfull \n");
        }
      for(int i=0; i<total_images;i++)
       { imagearray[i].initialize(ipm);
       }
   }
  else
   {  puts("\n New z-resolution cannot be zero");
      exit;
   }

  printf("\n Total no of images created are %d \n",total_images);
  return true;
}

/************************ inputarray****************************/
/* Input array is the struture which only stires the input images. The size of the array is equal to the no of input images*/

void allocateInputArray (int argc, char**argv, IPM & ipm)
{ inputarray = (Image*) malloc (sizeof(Image)*(argc-1));

  for (int i=0; i <argc-1; i++)

        { inputarray[i].initialize(ipm);
          inputarray[i].readImage (i+1, argv);
        }
}

/********************** paddedarray*****************************/
/* This array has a size of imagearray +2 in whhich the first and the last image has all the pixel defined as null pixel. Marching cube is only done on paddedarray. Here lenght denotes the original lenght of the array to be transfered to paddedarray */

// int  reallocate_paddedarray(Image array[], int length1, int length2)
//{if (length1>length2)
//  { for (int i = length2+2; i<length1; i++)
//      { delete[] paddedarray[i].imagedata;
//      }
//  }
// paddedarray=(Image*)realloc(paddedarray, sizeof(Image)*(length2+2));
// int columns = array[0].numberOfColumns;
// int rows = array[0].numberOfRows;
// int maxVal = array[0].maxVal;
// double x_z_pixel = array[0].x_pixelSize;
// double y_pixel = array[0].y_pixelSize;
// for (int i =0; i<(length2+2); i++)
//    { paddedarray[i].numberOfColumns =columns;
//      paddedarray[i].numberOfRows =rows;
//      paddedarray[i].maxVal= maxVal;
//      paddedarray[i].x_pixelSize =x_z_pixel;
//      paddedarray[i].z_pixelSize =x_z_pixel;
//      paddedarray[i].y_pixelSize =y_pixel;
//      paddedarray[i].reallocate(columns, rows);
//      if(i!=0 && i!=(length2+1))
//      { paddedarray[i].imagedata = array[i-1].imagedata;
//        paddedarray[i].paddImageBorders(); 
//      }
//      if (i==0|| i==(length2+1))
//      { for(long j =0; j< paddedarray[0].stringsize; j++)
//           {
//            paddedarray[i].imagedata[j]=0;
//           }
//      }
//    }
//  cout<<"\n exiting reallocation of paddedarray\n";
//  return (length2+2);
//}

/*********************Free all arrays **************************/

void  deconstruct_arrays (int argc, int length)
{  /*  free struct */
    for (int i=0; i<total_images; i++)
        { delete[] imagearray[i].imagedata;
        }
    for(int i=0; i<argc-1; i++)
        { delete[] inputarray[i].imagedata;
        }
    free(imagearray);
	if (imagearray = NULL)
	{
		puts("\n Imagearray  deallocation successfull \n"; );
	};
    free( inputarray);
	if (inputarray = NULL)
	{
		puts("\n Inputarray  deallocation successfull \n";);
	};
    
}


/****************************************************************/
/************* Image arrays processing functions****************/
/***************************************************************/

/* Interpolate used for pixelwise linear interpolation of images*/

void interpolate ( double & y_scan, double & y_AM, int & argc )
{ puts("\n Entered interpolation function \n");
  for( long i=0; i<(imagearray[0].getStringSize());i++)
  {  unsigned char  pixeldifference =0;//difference in pixel value
     unsigned char basepixel =0;
     bool pixelpositive = true;
     int inputno = 0; // tracking position in inputarray
     double yabsolute = 0.0; // tracking spatial position of the image
     double yreference = y_scan; // reference used to keep a check on interpolation
     double ycorrection =0.0; // correction to yreference after it moves to next consecutive pair of images
     int inpcheck = 1.0;  // counter used to avoid calculating slope in every iteration

    for(int j=0; j<total_images-2; j++)
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

	 if (yabsolute<= yreference && yreference <= (total_images-3)*y_AM)
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
	   if (yabsolute >yreference && yreference <= (total_images-3)*y_AM)
              { ycorrection = yreference;
		        yreference = yreference +y_scan;
	            inputno++;
                
	      }
      }
   } puts("\n Exiting Linear interpolation function ");
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
  for(long l =0; l<imagearray[0].getStringSize();l++)
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
    for(int n =1; n<total_images-1; n++)
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







/***************************************************************/
/**************** Image processing functions *******************/
/***************************************************************/

//void smoothing(Image& img, double sigma,int kernalsize)
//{ Image test;
//  int columns=img.numberOfColumns;
//  int rows = img.numberOfRows;
//  int maxval = img.maxVal;
//  int xgint;
//  test.initialize(columns, rows, maxval);
//  double r;
//  double  s=2.0 * sigma * sigma;
//  double sum =0.0;
//  double xg =0.0;
//  int width = test.numberOfColumns;
//  double gKernal[kernalsize][kernalsize];
//  double matrix[kernalsize][kernalsize];
//  int k = kernalsize/2;
//  Ip2 point,matrixmaker;
//  long int  matrixindex;
//
// for(int x =-k; x<=k; x++)
//  {  for(int z= -k; z<= k; z++)
//     {  r= sqrt(x*x+z*z);
//        gKernal[x+k][z+k]= (exp(-(r*r)/s))/(M_PI*s);
//        sum+=gKernal[x+k][z+k];
//     }
//  }
// cout<<"\n Kernal sum is "<<sum<<"\n";
// for(int i =0; i<kernalsize; i++)
//   { for(int j=0; j <kernalsize; j ++)
//     {  gKernal[i][j]/=sum;
//        cout<< gKernal[i][j]<<"\t";;
//     } cout<<"\n";
//   }
//
//
// for (long int i =0; i< test.stringsize-1; i++)
//    {xg=0; 
//     point = img.stringToCords(i);
//     cout<<(int) img.imagedata[i];
//     cout<< "\n";
//     if(point.x <rows-k && point.z<columns-k && point.x>k && point.z>k)
//       {   for(int x =-k; x<=k; x++)
//             {  for(int z= -k; z<= k; z++)
//                    { matrixmaker.x = point.x +x;
//                      matrixmaker.z = point.z +z;
//                      matrixindex =img.cordsToString(matrixmaker);
//                      matrix[x][z] =img.imagedata[matrixindex];
//                      cout<< matrix[x][z]<<"\t";
//		      xg += gKernal[x][z]*matrix[x][z];
//                    }
//              }
//       /*   for(int x =-k; x<=k; x++)
//              {  for(int z= -k; z<= k; z++)
//                    {  xg += gKernal[x][z]*matrix[x][z];
//                    }
//              } */
//
//        }
//        xgint = xg;
//        test.imagedata[i] =(unsigned char) xgint;
//        cout<<"\t\t"<<xg <<"\t"<<xgint <<"\n";
//        xg =0;
//    }
//
//  for (long int i=0; i< test.stringsize; i++)
//     { img.imagedata[i] = test.imagedata[i];
//     }
// }
//
//
//
////Threshold on original image
//
//void threshold( Image& img , unsigned char tLevel)
//{ 
//   for (long int j =0; j< img.stringsize; j++)
//     {  if ((tLevel+(unsigned char)5) <= img.imagedata[j] || (tLevel-(unsigned char)5) <= img.imagedata[j] )
//	  { img.imagedata[j] =(unsigned char) 255;
//	  }
//        else
//	   { img.imagedata[j]=(unsigned char) 0;
//	   }
//     }
//}

// Threshold and generate new image
void thresholdseries (unsigned char& tLevel)
{
	for (int i = 0; i < total_images; i++)
	{
		for (long int j = 0; j < imagearray[i].getStringSize(); j++)
		{
			if (tLevel <= imagearray[i].imagedata[j])
			{
				imagearray[i].imagedata[j] = (unsigned char)255;
			}
			else
			{
				imagearray[i].imagedata[j] = (unsigned char)0;
			}
		}
		
	}
}

// Cropping the image

void cropImage (Image& img,int & x_offset, int & z_offset, IPM &ipm)
{ Image temp;
  temp.initialize(ipm);
  Ip2 check, doublecheck;
  Ip2 startpoint, endpoint;
  char outimgname[20];
  long counter =0;
  long startindex;
  long endindex;
  
  startpoint.x = x_offset;
  startpoint.z = z_offset;
  startindex = img.cordsToString(startpoint);
  endpoint.x = startpoint.x + ipm.columns - 1;
  endpoint.z = startpoint.z + ipm.rows - 1;
  endindex = img.cordsToString(endpoint);

	  for (long i = startindex; i <= endindex; i++)
	  {
		  check = img.stringToCords(i);
		  if (check.x >= startpoint.x && check.x <= endpoint.x && check.z >= startpoint.z  && check.z <= endpoint.z)
		  {
			  temp.imagedata[counter] = img.imagedata[i];
			  counter++;
			  doublecheck.x = check.x; doublecheck.z = check.z;
		  }
	  }

	  img.reallocate(ipm.columns, ipm.rows);
	  for (long i = 0; i < img.getStringSize(); i++)
	  {
		  img.imagedata[i] = temp.imagedata[i];
	  }
    puts( "\n Done cropping the images ");
}

// Invert the image
void inversion(Image& temp)
{ unsigned char val; 
  for (long int j =0; j<temp.getStringSize(); j++)
     {  val = temp.imagedata[j];
        temp.imagedata[j]= temp.getMaxVal() -val;
     }
}





/***************************************************************/
/*********************** main **********************************/
/***************************************************************/

 int main( int argc, char** argv)
{  bool structarray;
//   unsigned char  Tlevel=0;
//   int t,paddedarraylength, kernalsize;
//   short crop,smooth, invert;
//   double y_scan, y_AM, x_z_pixel, sigma;
//   char outimgname[256];
//   long nTriang;
//   int columns, rows, maxVal,xCropOffset, zCropOffset;
//   cout <<"\n Enter the current pixel size in mm \n";
//   cin>>x_z_pixel;
//   cout<<"\n Enter the current y-resolution in mm \n";
//   cin>>y_scan;
//   cout<<"\n Enter the desired y-resolution in mm\n";
//   cin>>y_AM;
//   cout<<"\n Enter Threshold value \n";
//   cin >> t;
//   Tlevel = (unsigned char)t;
//   cout <<(int) Tlevel;
//   preInitializeCheck(columns, rows, maxVal, argv);
//
//   allocateInputArray ( argc, argv, columns, rows, maxVal, x_z_pixel,y_scan);
//
//
///*   cout<<"\n Do wou want to smooth the images ? \n Enter 1 for yes \n 2 for No \n";
//   cin>>smooth;
//   if (smooth ==1)
//   {  cout <<"\n Enter kernal size (Must be an odd number) \n";
//      cin>> kernalsize;
//      cout<<"\n Enter sigma value \n";
//      cin>> sigma;
//      for(int i=0; i<argc-1; i++)
//          {smoothing(inputarray[i],sigma, kernalsize);
//          }
//   } */
//
//   for(int i=0; i<argc-1;i++)
//      { snprintf(outimgname, sizeof(outimgname), "output/input_%d.pgm", i+1001);
//        inputarray[i].writeImage(outimgname);
//      }
//    
//   cout<<" \n Do you want to crop image? \n Enter \n 1 for YES \n 2 for NO \n";
//   cin>>crop;
//   if( crop == 1)
//   { cout<< "\n Enter the desired number of rows \n" ;
//     cin>> rows;
//     cout<<"\n Enter the desired number of columns \n";
//     cin >> columns;
//     cout<<"\n Enter the offset in x and z direction respectively \n";
//    cin >> xCropOffset >> zCropOffset;
//    for(int i=0; i<argc-1; i++)
//    {imageCrop(inputarray[i],rows, columns, xCropOffset, zCropOffset);
//    }
//   }
//
//   cout<<"\n Do you want to invert the colors of the image \n  Enter \n 1 for YES \n 2 for NO \n";
//   cin>> invert;
//   if (invert==1)
//    { for(int i=0; i<argc-1;i++)
//      { inversion(inputarray[i]);
//      }
//    }
// 
//   structarray = initialize_imagearray( y_scan, y_AM, argc , columns, rows, maxVal, x_z_pixel);
//   if(!structarray){ cout<<" \n Structure array was not initailized properly \n "; }
//  else 
//  { /*Switch between linear and Cubic by commenting one and uncommentig the other */
//     interpolate( y_scan, y_AM, argc);
//  // cubicInterpolate( y_scan, y_AM, argc); 
//
//  for(int i=0; i<total_images;i++)
//      { snprintf(outimgname, sizeof(outimgname), "output/Output_%d.pgm", i+1001);
//        imagearray[i].writeImage(outimgname);
//      }
//
///*  //threshold switch 
//    for(int i=0; i<total_images;i++)
//      { threshold(imagearray[i], Tlevel);
//      }
//*/
//
// // Stl of imagearray (output)
//  paddedarraylength =initialize_paddedarray(imagearray,total_images);
//  for(int i=0; i<paddedarraylength;i++)
//      {threshold(paddedarray[i] , Tlevel);
//        snprintf(outimgname, sizeof(outimgname), "output/padded1_%d.pgm", i+1001);
//        paddedarray[i].writeImage(outimgname);
//     }
//
//  nTriang = marchingCube (paddedarraylength,Tlevel);
//  snprintf(outimgname, sizeof(outimgname), "output/outstlbin.stl");
//  writeSTLbin(outimgname, nTriang);
//  cout <<" #triangles in output = "<< nTriang <<"\n ";
//
//  deconstruct_arrays(argc, paddedarraylength);
//  std::cout<<"\n Operation completed successfull ";
  return 0;
//}
}
