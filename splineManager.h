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
#include "shared.h"

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
  int nMode;
  float sysPar[NSYSMAX];
  float sysUnc[NSYSMAX];
  TTree* mcTree;
  fQreader* fqReader;

  //histograms
  TH1F* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TH1F* hMCNeut[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NMODE];
  //splines
  TSpine3* hSpline[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TSpline3 *hSplineNeut[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NMODE];
  //METHODS
  void setMCTree(TTree* tr);
  void setSysPar(int ipar, float val);
  void setSysUnc(int ipar, float val);
  float getEventWeight();
  void debugtest();
};

