#include "TH2Poly.h"

using namespace std;

////////////////////////////////////////////////////
// Class to make FV 2D histograms
//    X axis is towall
//    Y axis is wall
class TH2FV:public TH2Poly{

  public:

  TH2FV(const char* name, int bintype=0, int nbinsx =10, double xmin=0, double xmax=5000,
                                       int nbinsy =10, double ymin=0, double ymax=2000.);


  //vars
  int fBinType;
  int fNBinsX;
  int fNBinsY;
  double fXMin;
  double fXMax;
  double fYMin;
  double fYMax;

  //initialize
  void Init();

  protected:

  //for triangle plots:
  int LineIntersects(double* binx, double* biny);
  void SetSplitBins(double* splitx, double* splity);
  void InitSplitDiagonal();
  void InitFVBins();
  void InitFVBins2(); 
  void InitFVBins3();
};


//////////////////////////////////////////////////////
// FV bins for fits
void TH2FV::InitFVBins2(){

 //bin0
 double x0[] = {0,200,200,0};
 double y0[] = {0,0,200,0};
 AddBin(4,x0,y0);

 //bin1
// double x1[] = {200,400,400,200};
// double y1[] = {0,0,80,80};
// AddBin(4,x1,y1);

 //bin2
 double x2[] = {200,1000,1000,200};
 double y2[] = {0,0,80,80};
 AddBin(4,x2,y2);

 //bin3
 double x3[] = {1000,6000,6000,1000};
 double y3[] = {0,0,80,80};
 AddBin(4,x3,y3);

 //bin4
 double x4[] = {200,400,400,200};
 double y4[] = {80,80,200,200};
 AddBin(4,x4,y4);

 //bin5
 double x5[] = {400,1000,1000,400};
 double y5[] = {80,80,200,200};
 AddBin(4,x5,y5);

 //bin6
 double x6[] = {1000,6000,6000,1000};
 double y6[] = {80,80,200,200};
 AddBin(4,x6,y6);

 //bin7
 double x7[] = {200,6000,6000,200};
 double y7[] = {200,200,6000,200};
 AddBin(4,x7,y7);


 return;

}




//////////////////////////////////////////////////////
// FV bins for fits
void TH2FV::InitFVBins(){

 //bin0
 double x0[] = {0,200,200,0};
 double y0[] = {0,0,200,0};
 AddBin(4,x0,y0);

 //bin1
 double x1[] = {200,600,600,200};
 double y1[] = {0,0,80,80};
 AddBin(4,x1,y1);

 //bin2
 double x2[] = {600,6000,6000,600};
 double y2[] = {0,0,80,80};
 AddBin(4,x2,y2);

 //bin3
 double x3[] = {200,600,600,200};
 double y3[] = {80,80,200,200};
 AddBin(4,x3,y3);

 //bin4
 double x4[] = {600,6000,6000,600};
 double y4[] = {80,80,200,200};
 AddBin(4,x4,y4);

 //bin5
 double x5[] = {200,6000,6000,200};
 double y5[] = {200,200,6000,200};
 AddBin(4,x5,y5);

 return;

}



//////////////////////////////////////////////////////
// FV bins for perimeter
void TH2FV::InitFVBins3(){

 //bin0
 double x0[] = {0,300,300,0};
 double y0[] = {0,0,300,0};
 AddBin(4,x0,y0);

 //bin1
 double x1[] = {300,8000,8000,300};
 double y1[] = {0,0,80,80};
 AddBin(4,x1,y1);

 //bin2
 double x2[] = {8000,20000,20000,8000};
 double y2[] = {0,0,80,80};
 AddBin(4,x2,y2);

 //bin3
 double x3[] = {300,8000,8000,300};
 double y3[] = {80,80,300,300};
 AddBin(4,x3,y3);

 //bin4
 double x4[] = {8000,20000,20000,8000};
 double y4[] = {80,80,300,300};
 AddBin(4,x4,y4);

 //bin5
 double x5[] = {300,20000,20000,300};
 double y5[] = {200,200,6000,200};
 AddBin(4,x5,y5);

 return;

}



//////////////////////////////////////////////////////
// Initialize histogram with split diagonal bins
void TH2FV::InitSplitDiagonal(){

   //set some variables
   const int Nx = fNBinsX;
   const int Ny = fNBinsY;
   double xx[Nx+1];
   double yy[Ny+1];
   double wall = fYMin;
   double dtowall = (fXMax-fXMin)/(double)fNBinsX;
   double dwall = (fYMax-fYMin)/(double)fNBinsY;
   double towall = fXMin;

  double binxval[5];
  double binyval[5];


   //get coordinate arrays
   for (int i=0;i<=fNBinsX;i++){
     xx[i]=towall;
     towall+=dtowall;
   }
   for (int i=0; i<=fNBinsY; i++){
      yy[i] = wall;
      wall += dwall;
   }

   //loop over values and make bins
   for (int ix=0;ix<fNBinsX;ix++){
     for (int iy=0;iy<fNBinsY;iy++){
       binxval[0] = xx[ix];
       binyval[0] = yy[iy];
       binxval[1] = xx[ix+1];
       binyval[1] = yy[iy];
       binxval[2] = xx[ix+1];
       binyval[2] = yy[iy+1];
       binxval[3] = xx[ix];
       binyval[3] = yy[iy+1];
       
       //check to see if y=x intersects
       if (LineIntersects(binxval,binyval)){
         SetSplitBins(binxval, binyval); 
       }
       else{
         AddBin(4,binxval,binyval);
       }
     }
   }

   //return histogram
   return;
}


//////////////////////////////////////////////////////
// Adds split bins to histogram
void TH2FV::SetSplitBins(double* binxval, double* binyval){

  //is lower corner above or below line?
  int crosstype = 1; //< 1 for upper bin crossing
  if (binxval[0]>binyval[0]) crosstype = 2; //< 2 for lower bin crossing
  double splitbinx[5];
  double splitbiny[5]; 
  if (crosstype==1){
    //upper triagnle
    splitbinx[0]=binxval[0];
    splitbiny[0]=binxval[0];
    splitbinx[1]=binyval[2];
    splitbiny[1]=binyval[2];
    splitbinx[2]=binxval[3];
    splitbiny[2]=binyval[3];
    AddBin(3,splitbinx,splitbiny);
    //lower pentagon
    splitbinx[0]=binxval[0];
    splitbiny[0]=binyval[0];
    splitbinx[1]=binxval[1];
    splitbiny[1]=binyval[1];
    splitbinx[2]=binxval[2];
    splitbiny[2]=binyval[2];
    splitbinx[3]=binyval[2];
    splitbiny[3]=binyval[2];
    splitbinx[4]=binxval[3];
    splitbiny[4]=binxval[3];
    AddBin(5,splitbinx,splitbiny);
  }
  else{
    //lower triangle
    splitbinx[0]=binyval[0];
    splitbiny[0]=binyval[0];
    splitbinx[1]=binxval[1];
    splitbiny[1]=binyval[1];
    splitbinx[2]=binxval[1];
    splitbiny[2]=binxval[1];
    AddBin(3,splitbinx,splitbiny);
    //upper pentagon
    splitbinx[0]=binxval[0];
    splitbiny[0]=binyval[0];
    splitbinx[1]=binxval[0];
    splitbiny[1]=binyval[0];
    splitbinx[2]=binxval[1];
    splitbiny[2]=binxval[1];
    splitbinx[3]=binxval[2];
    splitbiny[3]=binyval[2];
    splitbinx[4]=binxval[3];
    splitbiny[4]=binxval[3];
    AddBin(5,splitbinx,splitbiny);
  }

  return;

}

//////////////////////////////////////////////////////
// See if diagonal line intersects this bin
int TH2FV::LineIntersects(double* binx, double* biny){
  
  //check lowest corner
  int intersect = 1;
  if (binx[1]>biny[1]) intersect*=-1;
  if (binx[2]>biny[2]) intersect*=-1;
  if (intersect<0) return 1;
  return 0;

}


void TH2FV::Init(){

  if (fBinType==1) InitSplitDiagonal();
  if (fBinType==0) InitFVBins();
  if (fBinType==2) InitFVBins2();
  if (fBinType==3) InitFVBins3(); //< template for perimeter bins
  return;
}


TH2FV::TH2FV(const char* name, int bintype, int nbinsx , double xmin, double xmax,
                                       int nbinsy , double ymin, double ymax){

  TH2Poly();

  fBinType = bintype;
  fXMin = xmin;
  fXMax = xmax;
  fYMin = ymin;
  fYMax = ymax;
  fNBinsX = nbinsx;
  fNBinsY = nbinsy;

  SetName(name);

  Init(); 
}


