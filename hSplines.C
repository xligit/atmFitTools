#ifndef HSPLINES_C
#define HSPLINES_C
#include "TH1F.h"
#include "TH2F.h"
#include "TSpline.h"
#include "TString.h"
#include <iostream>
using namespace std;

class hSplines{
  public:
  //a class for managing splines for a histogram
  hSplines(TH1F* h, int nsyst, const char* name="");
  hSplines(const char* name=""){}
 
  //spline creation 
  void buildSpline(int ibin, int isyst,double* X, double*Y, int N);

  //set spline
  void setSpline(int ibin, int isyst, TSpline3 *spline);
  //create modified histogram
  TH1F* buildModHistoAllPar(int npars, float* systPars);
  void buildModHisto(int npar, float parval);

  //getters
  TH1F* getBaseHisto() {return baseHisto;}
  TH1F* getModHisto()  {return modHisto;}
  
  //drawing
  void draw2D(int npts, int isyst);
  void drawSpline(int ibin, int isyst); 
  //debug test
  void debugTest();

  int checkSum;

  TH1F* baseHisto;
  TH1F* modHisto;
  TH2F* drawHisto;
  TSpline3 *theSpline[500][10]; //number of bins* number of sys pars
  TString nameTag;
  int nSyst; //number of systematic pars
  int nHistoBins; //number of bins in histogram
};

void  hSplines::drawSpline(int ibin, int isyst){
  if (theSpline[ibin][0]==NULL){
    cout<<"no such spine exists!"<<endl;
    return;
  }
  theSpline[ibin][0]->Draw();
}

void  hSplines::draw2D(int npts,int isyst){
  cout<<"getting par max and min"<<endl;
  float parmin= theSpline[1][isyst]->GetXmin();
  float parmax= theSpline[1][isyst]->GetXmax(); 
  float xx;
  float parval;
  cout<<"delete prev histo"<<endl;
  if (drawHisto!=NULL) drawHisto->Delete();
  cout<<"make new histo histo"<<endl;
  drawHisto= new TH2F("hdraw","hdraw",nHistoBins,baseHisto->GetBinLowEdge(1),
                      (baseHisto->GetBinLowEdge(baseHisto->GetNbinsX()) + baseHisto->GetBinWidth(baseHisto->GetNbinsX())),
                      npts,parmin,parmax);
  cout<<"fill xy histogram"<<endl;
  cout<<"nbins"<<nHistoBins<<endl;
  for (int xbin=1;xbin<=drawHisto->GetNbinsX();xbin++){
    for (int ybin=1;ybin<drawHisto->GetNbinsY();ybin++){
      xx = baseHisto->GetXaxis()->GetBinCenter(xbin); 
      parval = drawHisto->GetYaxis()->GetBinCenter(ybin);
      cout<<"parameter value: "<<parval<<endl;
      cout<<"spline value: "<< theSpline[xbin-1][isyst]->Eval(parval)<<endl;
      drawHisto->SetBinContent(xbin,ybin,theSpline[xbin-1][isyst]->Eval(parval)*baseHisto->GetBinContent(xbin)); 
    }
  }
  cout<<"draw histogram"<<endl;
  drawHisto->Draw("colz");
  return;
}

void  hSplines::debugTest(){
  //test it
  double x1[5] = {0.,1.,2.,3.,4.};
  double y1[5] = {0.5,0.1,0.9,1.3,1.1};
  nSyst = 1;
  nHistoBins = 5;
  baseHisto = new TH1F("hdebug","hdebug",nHistoBins,0,5);
  baseHisto->Fill(1.);
  baseHisto->Fill(1.);
  baseHisto->Fill(2.);
  baseHisto->Fill(3.);
  baseHisto->Fill(4.);
  baseHisto->Fill(5.);
  baseHisto->Fill(0.);
  baseHisto->Fill(4.);
  baseHisto->Fill(4.);
  baseHisto->Fill(4.);
  for (int ibin=1;ibin<=nHistoBins;ibin++){
    buildSpline(ibin,0,x1,y1,5);
  }
  draw2D(10,0);
  return;
}

void hSplines::buildModHisto(int ipar, float parval){
  if (modHisto==NULL){   
    TString modhname = baseHisto->GetName();
    cout<<"creating mod histo from base: "<<baseHisto->GetName()<<endl;
    modhname.Append("_modified");
//    int basenbins = baseHisto->GetNbinsX();
  //  float basexmin  = baseHisto->GetBinLowEdge(1);
   // float  basexmax = baseHisto->GetBinLowEdge(basenbins);
    modHisto=(TH1F*)baseHisto->Clone(modhname.Data());
  }
  float newcontent;
  for (int ibin=1;ibin<=nHistoBins;ibin++){
    newcontent = baseHisto->GetBinContent(ibin);
    if (theSpline[ibin][ipar]->Eval(parval)==0){
      cout<<"zero spline interpolation for bin"<<ibin<<" with content "<<newcontent<<endl;
    }
    newcontent *= theSpline[ibin][ipar]->Eval(parval); 
    modHisto->SetBinContent(ibin,newcontent);
  }
  return;
}

TH1F* hSplines::buildModHistoAllPar(int npars, float *systPars){
  float binx;
  float newcontent;
 if (modHisto==NULL){   
    TString modhname = baseHisto->GetName();
    cout<<"creating mod histo from base: "<<baseHisto->GetName()<<endl;
    modhname.Append("_modified");
    modHisto=(TH1F*)baseHisto->Clone(modhname.Data());
  }
  cout<<"nbins: "<<nHistoBins<<endl;
  for (int ibin=1;ibin<=nHistoBins;ibin++){
    newcontent = baseHisto->GetBinContent(ibin);
    binx = baseHisto->GetBinCenter(ibin);
    for (int isyst=0;isyst<npars;isyst++){
      newcontent *= theSpline[ibin][isyst]->Eval(systPars[isyst]); 
    }
    modHisto->SetBinContent(ibin,newcontent);
  }
  return modHisto;
}

void hSplines::setSpline(int ibin, int isyst, TSpline3 *spline){
  theSpline[ibin][isyst] = spline;
  checkSum--;
  return;
}

void hSplines::buildSpline(int ibin, int isyst,double* X, double*Y, int N){
 TString splineName = nameTag.Data();
 splineName.Append(Form("_spline_bin%d_par%d",ibin,isyst));
 theSpline[ibin][isyst] = new TSpline3(splineName.Data(),X,Y,N); 
 checkSum--;
 return;
}

hSplines::hSplines(TH1F* h, int nsyst, const char* name){
  if (name!=NULL) nameTag = name;
   nameTag = h->GetName();
   nHistoBins=h->GetNbinsX();
   nSyst=nsyst;
   baseHisto=h;
   TString modhname = h->GetName();
   cout<<"base histogram name: "<<h->GetName()<<endl;
   modhname.Append("_modified");
   int basenbins = baseHisto->GetNbinsX();
   float basexmin  = baseHisto->GetBinLowEdge(1);
   float  basexmax = baseHisto->GetBinLowEdge(basenbins);
   checkSum = nHistoBins*nSyst;

   //modHisto= new TH1F(modhname.Data(),modhname.Data(),basenbins,basexmin,basexmax);
   
//   modHisto=(TH1F*)baseHisto->Clone(modhname.Data());
//   modHisto->Reset();
}

#endif
