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
  TH1F* buildModHisto(int npars, float* systPars);

  //getters
  TH1F* getBaseHisto() {return baseHisto;}
  TH1F* getModHisto()  {return modHisto;}
  
  //drawing
  void draw2D(int npts, int isyst);
  void drawSpline(int ibin, int isyst); 
  //debug test
  void debugTest();

  private:
  TH1F* baseHisto;
  TH1F* modHisto;
  TH2F* drawHisto;
  TSpline3 *theSpline[500][10]; //number of bins* number of sys pars
  TString nameTag;
  int nSyst; //number of systematic pars
  int nBins; //number of bins in histogram
};

void  hSplines::drawSpline(int ibin, int isyst){
  if (theSpline[ibin][0]==NULL){
    cout<<"no such spine exists!"<<endl;
    return;
  }
  theSpline[ibin][0]->Draw();
}

void  hSplines::draw2D(int npts,int isyst){
  float parmin= theSpline[1][isyst]->GetXmin();
  float parmax= theSpline[1][isyst]->GetXmax(); 
  float xx;
  float parval;
  if (drawHisto) drawHisto->Delete();
  drawHisto= new TH2F("hdraw","hdraw",nBins,baseHisto->GetBinLowEdge(1),
                      (baseHisto->GetBinLowEdge(baseHisto->GetNbinsX()) + baseHisto->GetBinWidth(baseHisto->GetNbinsX())),
                      npts,parmin,parmax);
  for (int xbin=1;xbin<=drawHisto->GetNbinsX();xbin++){
    for (int ybin=1;ybin<drawHisto->GetNbinsY();ybin++){
      xx = baseHisto->GetXaxis()->GetBinCenter(xbin); 
      parval = drawHisto->GetYaxis()->GetBinCenter(ybin);
      drawHisto->SetBinContent(xbin,ybin,(baseHisto->GetBinContent(xbin)*theSpline[xbin][isyst]->Eval(parval))); 
    }
  }
  drawHisto->Draw("colz");
  return;
}

void  hSplines::debugTest(){
  //test it
  double x1[5] = {0.,1.,2.,3.,4.};
  double y1[5] = {0.5,0.1,0.9,1.3,1.1};
  nSyst = 1;
  nBins = 5;
  baseHisto = new TH1F("hdebug","hdebug",nBins,0,5);
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
  for (int ibin=1;ibin<=nBins;ibin++){
    buildSpline(ibin,0,x1,y1,5);
  }
  draw2D(10,0);
  return;
}

TH1F* hSplines::buildModHisto(int npars, float *systPars){
  float binx;
  float newcontent;
  for (int ibin=1;ibin<=nBins;ibin++){
    newcontent = baseHisto->GetBinContent(ibin);
    binx = baseHisto->GetBinCenter(ibin);
    for (int isyst=0;isyst<nSyst;isyst++){
      newcontent *= theSpline[ibin][isyst]->Eval(binx); 
    }
    modHisto->SetBinContent(ibin,newcontent);
  }
  return modHisto;
}

void hSplines::setSpline(int ibin, int isyst, TSpline3 *spline){
  theSpline[ibin][isyst] = spline;
  return;
}

void hSplines::buildSpline(int ibin, int isyst,double* X, double*Y, int N){
 TString splineName = nameTag.Data();
 splineName.Append(Form("_spline_bin%d_par%d",ibin,isyst));
 cout<<"creating spline: "<<splineName.Data()<<endl;
 theSpline[ibin][isyst] = new TSpline3(splineName.Data(),X,Y,N); 
 return;
}

hSplines::hSplines(TH1F* h, int nsyst, const char* name){
  if (name!=NULL) nameTag = name;
   nameTag = h->GetName();
   nBins=h->GetNbinsX();
   nSyst=nsyst;
}
