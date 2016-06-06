#include "TH2Poly.h"
#include "TVector2.h"
#include "TMath.h"


#define NHBINSMAX 9999

using namespace std;

////////////////////////////////////////////////////
// Class to make FV 2D histograms
//    X axis is towall
//    Y axis is wall
class TH2FV:public TH2Poly{

  public:

  TH2FV(const char* name, int bintype=0, int nbinsx =10, double xmin=0, double xmax=5000,
                                       int nbinsy =10, double ymin=0, double ymax=2000.);
//  TH2FV(const char* name, int bintype, int nbinsx, double xmin, double xmax,
                                       //int nbinsy, double ymin, double ymax);





  //vars
  int fBinType;
  int fNBinsX;
  int fNBinsY;
  double fXMin;
  double fXMax;
  double fYMin;
  double fYMax;
  double fCenter[NHBINSMAX][2];

  //initialize
  void Init();
  
  //add bin and keep track of center
  //(messes up AddBin() for some reason
  void AddBinWithCenter(int n, double *xx, double* yy);

  //add bin center
  void AddBinCenter(int n, double *xx, double* yy);

  // get the center of a bin
  double GetBinCenterX(int nbin);
  double GetBinCenterX(int binx, int biny);
  double GetBinCenterY(int nbin);
  double GetBinCenterY(int binx, int biny);
  double GetMaxWall();

  // draw standard view
  void DrawStdView(const char* opts);
 
  // test bin consistency
  void TestHisto();

  double fMaxWall;

  protected:

  //for triangle plots:
  int LineIntersects(double* binx, double* biny);
  //double fMaxWall;
  void SetSplitBins(double* splitx, double* splity);
  void InitSplitDiagonal();
  void InitFVBins();
  void InitFVBins2(); 
  void InitFVBins3();
  void InitFVBins4();
  void InitFVBins5();
  void InitFVBins6();
  void InitFVBins7();
  void InitStdBins(double wall1, double wall2, double towall1,
                   double towall2, double towall3, double towall4);
  
};

double TH2FV::GetMaxWall(){
  return fMaxWall;
}

void TH2FV::InitStdBins(double wall1, double wall2, double towall1,
                   double towall2, double towall3, double towall4){




// double wall1 = 80.;
// double wall2 = 200.;
// double towall1 = 300.;
// double towall2 = 6000.;
// double towall3 = 800.;
// double towall4 = 800.;
 fMaxWall = towall4;

 //bin1
 double x0[] = {0,towall1,towall1,0};
 double y0[] = {0,0,towall1,0};
// AddBin(4,x0,y0);
// AddBinCenter(3,x0,y0);
 AddBinWithCenter(4,x0,y0);

 //bin2
 double x2[] = {towall1,towall2,towall2,towall1};
 double y2[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x2,y2);

 //bin3
// double x3[] = {towall2,6000,6000,towall2};
// double y3[] = {0,0,wall1,wall1};
// AddBinWithCenter(4,x3,y3);

 //bin4
 double corner1 = TMath::Max(towall1,wall2);
 double corner2 = TMath::Min(towall1,wall2);
 double x4[] = {towall1, towall3,towall3, corner1, towall1,towall1};
 double y4[] = {wall1,   wall1,   wall2,  wall2,  corner2, wall1};
 AddBinWithCenter(6,x4,y4);

 //bin5
 double x5[] = {towall3,6000,6000,towall3};
 double y5[] = {wall1,wall1,wall2,wall2};
 AddBinWithCenter(4,x5,y5);

 //bin6
 corner1 = TMath::Max(towall1,wall2);
 double x6[] = {corner1,towall4,towall4,corner1,corner1};
 double y6[] = {wall2,  wall2,  towall4,corner1,wall2};
 AddBinWithCenter(5,x6,y6);

 //bin7
 //double x7[] = {towall4,6000,6000,towall4};
 //double y7[] = {wall2,wall2,towall4,towall4};
 double x7[] = {towall4,6000,6000,towall4};
 double y7[] = {wall2, wall2,6000,towall4};
 AddBinWithCenter(4,x7,y7);

 //bin8
// double xx[] = {towall4, 6000.,6000.,towall4};
// double yy[] = {towall4,towall4,6000., towall4};
// double xx[]={1000,2000,2000,1000};
// double yy[]={2000,2000,3000,3000};
// AddBinWithCenter(4,xx,yy);
 


 return;


}

void TH2FV::TestHisto(){
  
  for (int i=1; i<=GetNumberOfBins(); i++){
    double xx = GetBinCenterX(i);
    double yy = GetBinCenterY(i);
    //cout<<"Fill: "<<"("<<xx<<","<<yy<<")"<<endl;
    Fill(xx,yy);
  }

  Draw("colz");
}

void TH2FV::DrawStdView(const char* opts){
  GetXaxis()->SetRangeUser(0,1600);
  GetYaxis()->SetRangeUser(0,1400);
  Draw(opts);
  return;
}

double TH2FV::GetBinCenterY(int binx, int biny){
  int globalbin = GetBin(binx, biny);
  return fCenter[globalbin-1][1];
}

double TH2FV::GetBinCenterY(int nbin){
  return fCenter[nbin-1][1];
}

double TH2FV::GetBinCenterX(int binx, int biny){
  int globalbin = GetBin(binx, biny);
  return fCenter[globalbin-1][0];
}

double TH2FV::GetBinCenterX(int nbin){
  return fCenter[nbin-1][0];
}


//////////////////////////////////////////////////////
// Adds a bin and updates the list of bin centers
void TH2FV::AddBinWithCenter(int n, double* xx, double* yy){

  //Add the bin
  AddBin(n,xx,yy);

  //first calculate the center
  TVector2 vcenter;
  vcenter.Set(0.,0.);
  for (int i=0; i<n; i++){
    TVector2 vertex;
    vertex.Set(xx[i],yy[i]);
    vcenter = vcenter + vertex;
  }
  double norm = 1./(double)n;
  vcenter = vcenter*norm;

  // fill the center array
  int N = GetNumberOfBins();
  fCenter[N-1][0] = vcenter.X();
  fCenter[N-1][1] = vcenter.Y();

  return;

}

//////////////////////////////////////////////////////
// just adds bin center to list
void TH2FV::AddBinCenter(int n, double* xx, double* yy){

  //first calculate the center
  TVector2 vcenter;
  vcenter.Set(0.,0.);
  for (int i=0; i<n; i++){
    TVector2 vertex;
    vertex.Set(xx[i],yy[i]);
    vcenter = vcenter + vertex;
  }
  double norm = 1./(double)n;
  vcenter = vcenter*norm;

  // fill the center array
  int N = GetNumberOfBins();
  fCenter[N-1][0] = vcenter.X();
  fCenter[N-1][1] = vcenter.Y();

  return;

}



//////////////////////////////////////////////////////
// FV bins for fits
void TH2FV::InitFVBins6(){


 double wall1 = 80.;
 double wall2 = 200.;
 double towall1 = 500.;
 double towall2 = 1100.;
 double towall3 = 1100.;
 double towall4 = 1100.;

 //bin0
 double x0[] = {0,towall1,towall1,0};
 double y0[] = {0,0,towall1,0};
 AddBinWithCenter(4,x0,y0);

 //bin1
 double x2[] = {towall1,towall2,towall2,towall1};
 double y2[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x2,y2);

 //bin2
 double x3[] = {towall2,6000,6000,towall2};
 double y3[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x3,y3);

 //bin3
 double x4[] = {towall1,towall3,towall3,towall1};
 double y4[] = {wall1,wall1,wall2,wall2};
 AddBinWithCenter(4,x4,y4);

 //bin4
 double x5[] = {towall3,6000,6000,towall3};
 double y5[] = {wall1,wall1,wall2,wall2};
 AddBinWithCenter(4,x5,y5);

 //bin5
 double x6[] = {towall1,towall4,towall4,towall1,towall1};
 double y6[] = {wall2,wall2,towall4,towall1,wall2};
 AddBinWithCenter(4,x6,y6);

 //bin6
 double x7[] = {towall4,6000,6000,towall4};
 double y7[] = {wall2,wall2,towall4,towall4};
 AddBinWithCenter(4,x7,y7);

 //bin7
 double x8[] = {towall4,6000,6000,towall4};
 double y8[] = {towall4,towall4,6000,towall4};
 AddBinWithCenter(4,x8,y8);
 


 return;

}



//////////////////////////////////////////////////////
// FV bins for fits
void TH2FV::InitFVBins7(){


 double wall1 = 140.;
 double wall2 = 260.;
 double towall1 = 500.;
 double towall2 = 1100.;
 double towall3 = 1100.;
 double towall4 = 1100.;




 //bin0
 double x0[] = {0,towall1,towall1,0};
 double y0[] = {0,0,towall1,0};
 AddBinWithCenter(4,x0,y0);

 //bin1
 double x2[] = {towall1,towall2,towall2,towall1};
 double y2[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x2,y2);

 //bin2
 double x3[] = {towall2,6000,6000,towall2};
 double y3[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x3,y3);

 //bin3
 double x4[] = {towall1,towall3,towall3,towall1};
 double y4[] = {wall1,wall1,wall2,wall2};
 AddBinWithCenter(4,x4,y4);

 //bin4
 double x5[] = {towall3,6000,6000,towall3};
 double y5[] = {wall1,wall1,wall2,wall2};
 AddBinWithCenter(4,x5,y5);

 //bin5
 double x6[] = {towall1,towall4,towall4,towall1,towall1};
 double y6[] = {wall2,wall2,towall4,towall1,wall2};
 AddBinWithCenter(4,x6,y6);

 //bin6
 double x7[] = {towall4,6000,6000,towall4};
 double y7[] = {wall2,wall2,towall4,towall4};
 AddBinWithCenter(4,x7,y7);

 //bin7
 double x8[] = {towall4,6000,6000,towall4};
 double y8[] = {towall4,towall4,6000,towall4};
 AddBinWithCenter(4,x8,y8);
 


 return;

}





//////////////////////////////////////////////////////
// FV bins for fits
void TH2FV::InitFVBins4(){


 double wall1 = 80.;
 double wall2 = 200.;
 double towall1 = 300.;
 double towall2 = 6000.;
 double towall3 = 800.;
 double towall4 = 800.;

 InitStdBins(wall1,wall2,towall1,towall2,towall3,towall4);

 return;

}



//////////////////////////////////////////////////////
// FV bins for fits
void TH2FV::InitFVBins5(){


 double wall1 = 140.;
 double wall2 = 260.;
 double towall1 = 200.;
 double towall2 = 800.;
 double towall3 = 800.;
 double towall4 = 800.;

 //bin0
 double x0[] = {0,towall1,towall1,0};
 double y0[] = {0,0,towall1,0};
 AddBinWithCenter(4,x0,y0);

 //bin1
 double x2[] = {towall1,towall2,towall2,towall1};
 double y2[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x2,y2);

 //bin2
 double x3[] = {towall2,6000,6000,towall2};
 double y3[] = {0,0,wall1,wall1};
 AddBinWithCenter(4,x3,y3);

 //bin3
 double x4[] = {towall1,towall3,towall3,wall2,towall1};
 double y4[] = {wall1,wall1,wall2,wall2,towall1};
 AddBinWithCenter(5,x4,y4);

 //bin4
 double x5[] = {towall3,6000,6000,towall3};
 double y5[] = {wall1,wall1,wall2,wall2};
 AddBinWithCenter(4,x5,y5);

 //bin5
 double x6[] = {wall2,towall4,towall4,wall2};
 double y6[] = {wall2,wall2,towall4,wall2};
 AddBinWithCenter(4,x6,y6);

 //bin6
 double x7[] = {towall4,6000,6000,towall4};
 double y7[] = {wall2,wall2,towall4,towall4};
 AddBinWithCenter(4,x7,y7);

 //bin7
 double x8[] = {towall4,6000,6000,towall4};
 double y8[] = {towall4,towall4,6000,towall4};
 AddBinWithCenter(4,x8,y8);
 


 return;

}



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

  double binxval[4];
  double binyval[4];


   //get coordinate arrays
   for (int i=0;i<=fNBinsX;i++){
     xx[i]=towall;
 //    cout<<"xx"<<i<<" = "<<towall<<endl;
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
         AddBinWithCenter(4,binxval,binyval);
        // AddBin(4,binxval,binyval);
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
    AddBinWithCenter(3,splitbinx,splitbiny);
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
    AddBinWithCenter(5,splitbinx,splitbiny);
  }
  else{
    //lower triangle
    splitbinx[0]=binyval[0];
    splitbiny[0]=binyval[0];
    splitbinx[1]=binxval[1];
    splitbiny[1]=binyval[1];
    splitbinx[2]=binxval[1];
    splitbiny[2]=binxval[1];
    AddBinWithCenter(3,splitbinx,splitbiny);
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
    AddBinWithCenter(5,splitbinx,splitbiny);
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

  if (fBinType==-1) InitSplitDiagonal();
  if (fBinType==0) InitFVBins();
  if (fBinType==1) InitStdBins(80.,200.,300.,6000.,800.,800.);
  if (fBinType==2) InitStdBins(80.,200.,550.,6000.,1050.,1050.);
  if (fBinType==3) InitStdBins(140.,260.,300.,6000.,800.,800.);
  if (fBinType==4) InitStdBins(140.,260.,550.,6000.,1050.,1050.);




//  if (fBinType==2) InitFVBins2();
//  if (fBinType==3) InitFVBins3(); //< template for perimeter bins
//  if (fBinType==4) InitFVBins4(); //< template for perimeter bins
//  if (fBinType==5) InitFVBins5(); //< template for perimeter bins
//  if (fBinType==6) InitFVBins6(); //< template for perimeter bins
//  if (fBinType==7) InitFVBins7(); //< template for perimeter bins



  return;
}


TH2FV::TH2FV(const char* name, int bintype, int nbinsx , double xmin, double xmax,
                                       int nbinsy , double ymin, double ymax){

  TH2Poly();
  GetXaxis()->SetTitle("Towall [cm]");
  GetYaxis()->SetTitle("Wall [cm]");
  GetYaxis()->SetTitleOffset(1.2);
  fBinType = bintype;
  fXMin = xmin;
  fXMax = xmax;
  fYMin = ymin;
  fYMax = ymax;
  fNBinsX = nbinsx;
  fNBinsY = nbinsy;

  SetName(name);

  ChangePartition(100,100);

  Init(); 
}


