#include"CTimage.h"

using namespace std; 
// In class functions 
//Default constructer 
//Image::Image()
//{ this->numberofcolumns =512;
//  this->numberofrows =512;
//  this->x_pixelsize =0.4;
//  this->z_pixelsize =0.4;
//  this->maxval =255;
//  this-> stringsize = ((this->numberofrows)*(this->numberofcolumns));
//  this->imagedata  =(unsigned char*) malloc (sizeof(unsigned char)*(stringsize));
//  this->imageheader  =( char*) malloc (sizeof(char)*(15));
//  if(this->imagedata==null)
//    { cout<<"\n  space allocation for imagedata was not sucessfull ";
//    }
//}


//Dynamic constructer (preferred)
void Image::initialize(const IPM& ipm)
{numberOfColumns =ipm.columns;
 numberOfRows =ipm.rows;
 maxVal=ipm.maxval;
 x_z_pixelSize =ipm.x_z_pixel;
 y_pixelSize =ipm.y_pixel;
 stringsize = ((this->numberOfRows)*(this->numberOfColumns));
 imagedata  =(unsigned char*) malloc (sizeof(unsigned char)*(stringsize));
 imageheader  =( char*) malloc (sizeof(char)*(15));
 if(imagedata==NULL)
    { puts("\n  space allocation for Imagedata was not sucessfull ");
    }
}

/* Input : Desired no of rows and columns to change the size of imagedata */
void Image::reallocate(const int& columns,const int& rows)
{numberOfColumns =columns;
 numberOfRows =rows;
 stringsize = ((numberOfRows)*(numberOfColumns));
 imagedata  =(unsigned char*) realloc (imagedata, sizeof(unsigned char)*(stringsize));
  if(imagedata==NULL)
    { puts("\n  space reallocation for Imagedata was not sucessfull ");
    }
}

/*************************************************************/
/****************** read image functions**********************/
/*************************************************************/

void Image::readHeader (FILE *fp)
{
	char magicNumber = '0';
	fscanf(fp, " P%c %d %d %d",&magicNumber, &numberOfColumns, &numberOfRows, &maxVal);
	if (magicNumber != '5')
	{ // replace with a suitable exception
		puts("\n Unsupported file ");
	}
	sprintf(this->imageheader, "P%c\n%d %d\n%d\n",&magicNumber, &numberOfColumns, &numberOfRows,&maxVal);
	if (this->imagedata == NULL) {
		fputs(" Memory error", stderr);
		exit(2);
	}
}

void  Image::readImage ( int fileno, char** argv) 
{
   long lSize; 
   int headerSize;
   size_t result; // to cross-check size of the imagedata
   int location =0;

   FILE *fp = fopen(argv[fileno],"rb");
//to replace with exceptions
   if(fp == NULL )
    { puts("\n Such a file does not exist \n");
      exit(1);;
    }
   else
    { puts("\n Successfully loaded the image");
    }
   // reads image header
   readHeader(fp);
   headerSize = sizeof(this->imageheader) - 1;
   fseek( fp,0,SEEK_END);
   lSize = ftell (fp)-headerSize;
   rewind(fp);
   fseek(fp,headerSize , SEEK_SET);
   result = fread(this->imagedata, sizeof(unsigned char), lSize , fp);
   if (result != lSize) {
	   fputs("Reading error", stderr); exit(3);
   }
   puts("\n Done reading the image.. \n");
   fclose(fp);
}

/***************************************************************/
/*****************write image functions**************************/
/****************************************************************/

// rewrites imageheader as size of the image might change
void  Image::writeHeader ()
{
  sprintf(this->imageheader, "P5\n%d %d\n%d\n",this->numberOfColumns,this->numberOfRows, this->maxVal);

}
 // Writes image to a file named by foutname[]
void  Image::writeImage (char foutname[])
{ FILE *fout;
 fout = fopen(foutname, "wb");
 if(fout == NULL )
    { puts("\n error while opening the fout \n");
      exit(1);;
    }

  writeHeader();
  fwrite(this->imageheader, sizeof(char),15, fout);
  fwrite("\n",sizeof(char),1,fout); 
  fwrite(this->imagedata, sizeof(unsigned char),this->stringsize, fout);
  if( fout ==NULL)
  { cout<<"\n Error writing the image \n ";
  }
  else
  {
  cout<<"\n Done writing the image \n";
  }
  fclose(fout);
}

/***************************************************************/
/******conversion functions between bitmap and imagedata*******/
/* Here the bitmap is assumed to start from 0 to N in both x and y directions    */
/*****  Also imagedata is an array so the index starts from 0 ***/ 
/****************************************************************/

   Ip2 Image:: stringToCords (const long int & i)
   { int x, z;
     x=(i % this->numberOfColumns);
     z= (i/this->numberOfColumns);
     return Ip2(x,z);
   }

  long  int Image::cordsToString (Ip2& p)
   { long int i;
     i = (p.x+p.z*this->numberOfColumns);
     return i;
   }


  Dp3 Image:: stringToRealCords (const long int& i, const int& yInt)
  { double x, y, z;
    x=(i % this->numberOfColumns)*this->x_z_pixelSize;
    z= (i/this->numberOfColumns)*this->x_z_pixelSize;
    y =(yInt * this->y_pixelSize);
    return Dp3(x,y,z);
  } 


/****************************************************************/
/************* additional functions *************************/
/*****************************************************************/
void Image:: paddImageBorders()
{ Ip2 check;
   for(long int i=0; i< this->stringsize-1; i++)
   { check = stringToCords (i);
     if (check.x == 0 || check.z==0 || check.z == (this-> numberOfRows-1) ||check.x == (this->numberOfColumns-1))
        { this-> imagedata[i]=0;
        }
   }
} 


//int realCordsToString (double, double, char*);
//   ~Image(){};
 