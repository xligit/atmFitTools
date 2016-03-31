#ifndef HSPLINES_H
#define HSPLINES_H

#include "TH1D.h"
#include "TH2F.h"
#include "TSpline.h"
#include "TString.h"
#include "shared.h"
#include <iostream>

using namespace std;

class hSplines{
  public:
  //a class for managing splines for a histogram
  hSplines(TH1D* h, int nsyst, const char* name="");

  hSplines(const char* name=""){}
  ~hSplines();
 
  //spline creation 
  void buildSpline(int ibin, int isyst,double* X, double*Y, int N);

  //set spline
  void setSpline(int ibin, int isyst, TSpline3 *spline);
  //create modified histogram
  TH1D* buildModHistoAllPar(int npars, double* systPars);
  void buildModHisto(int ipar, double parval);

  //getters
  TH1D* getBaseHisto() {return baseHisto;}
  TH1D* getModHisto()  {return modHisto;}
  
  //drawing
  void draw2D(int npts, int isyst);
  void drawSpline(int ibin, int isyst); 
  //debug test
  void debugTest();

  int checkSum;

  TH1D* baseHisto; //the default histogram from which all modified histograms are obtained
  TH1D* modHisto;  //pointer to the histogram that has been modified by splines
  TH2F* drawHisto;  //pointer to the histogram that is used in the draw functions
  TSpline3 *theSpline[NHBINSMAX][NSYSTMAX]; //number of bins* number of sys pars
  double evaluateSpline(int ibin, int ipar, double parvalue);
  TString nameTag;
  int nSyst; //number of systematic pars
  int nHistoBins; //number of bins in histogram
};

#endif
