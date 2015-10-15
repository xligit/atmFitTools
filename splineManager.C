#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "fQreader.C"
#include <iostream>
#include "THStack.h"
#include "TLegend.h"
#include "TSpline3.h"

#define NSAMPMAX 5
#define NCOMPMAX 20
#define NATTMAX 20
#define NBINMAX 10
#define NSYSMAX 10
//creates and manages splines for flux/xsec systematics
class splineManager{
  public:
  //CONSTRUCTOR
  splineManager(const char* name, int nsample, int nbin, int ncomp, int natt, int nsyst);

  //INTERNAL VARIABLES
  int nBin;
  int nComp;
  int nSamp;
  int nAtt;
  int nSyst;
  float sysPar[NSYSMAX];
  float sysUnc[NSYSMAX];
  TTree* mcTree;
  fQreader* fqReader;

  //histograms
  TH1F* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  //splines
  TSpine3* hSpline[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  //METHODS
  void setMCTree(TTree* tr);
  void setSysPar(int ipar, float val);
  void setSysUnc(int ipar, float val);
  float getEventWeight();
  void debugtest();
};

void getEventWeight(){
  float ww=0.;
  if (fqReader->mode==1)ww*=0.9;
  return ww;
}

void splineManager::splineManager(const char* name, int nsample, int nbin, int ncomp, int natt, int nsyst){
  int nBin=nbin;
  int nComp=ncomp;
  int nAtt=natt;
  int nSyst=nsyst;
}

void splineManager::setSysUnc(int ipar, float val){
  sysUnc[ipar]=val;
  return;
}

void splineManager::setSysPar(int ipar, float val){
  sysPar[ipar]=val
  return;
}

void splineManager::setMCTree(TTree* tr){
  mcTree = tr;
}






