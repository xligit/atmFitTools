#ifndef __T2KPREPROCESS_H__
#define __T2KPREPROCESS_H__

#include <iostream>
#include "visRing.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph.h"
#include "TString.h"
#include "FVCalculators.cxx"
#include "sharedPars.cxx"
#include "shared.h"
#include "t2kfqEvent.h"
#include "TH2FV.cxx"
#include "masktools.cxx"
#include "TObjArray.h"

using namespace std;

#define NFILEMAX 5000

//////////////////////////////////////////////////////////////
//Class to take care of preprocessing of all data and MC files
//Usage:
//  1) Specify parameters in sharedpars.dat file.  These parameters 
//     specify the files to be run on, where the output should go,
//     and which FV binning and MC component definitions to use
//  2) Call runPreprocess()

class t2kPreProcess{
  public:

  /////////////////////
  //constructors
  t2kPreProcess();
  t2kPreProcess(TChain* chin,const char* name="");
  t2kPreProcess(TTree* trin,const char* name="");
  t2kPreProcess(TChain *mc, TChain *spline, const std::string name);

  /////////////////////
  //internal variables
  TChain* chmc;
  TChain* chdat;
  TChain *chspline;
  TTree *trspline;
  TTree* tr;
  TTree* trout;
  TFile* fout;
  TString nameTag;
  TString outDir;
  TString fileNames[NFILEMAX];
  TString parFileName;
  t2kfqEvent* fq;
  visRing*  vis;
  int MCComponents;
  int FVBinning;
  int nFiles;
  int nFilesSpline;
  int MCSamples;

  TH1D *hWeight;
  TGraph *gWeight;
  std::map<std::string, double> attributeMap;
  std::string attributeList[50];
  int nAttributes;
  TH1D *hmask;
  bool isData;
  bool existSpline;
  ////////////////////////
  //for cuts
  int NHITACMax;
  double EVisMin;
  double WallMin;
  double ToWallMin;
  int NSEMax;
  int NSEMin;
  double InGateMin; 
  TH2FV* hFVBins;

  //////////////////////
  //new branch variables
  float towall;
  float wall;
  float evtweight; // banff weight
  float rfgweight;
  float attribute[1000];
  float towallv[50];
  float minconev[50];
  float perimv[50];
  float wallv2;
  float fq1rwall[10][7];
  float fq1rtowall[10][7];
  float fq1rmincone[10][7];
  float fq1rperim[10][7];
  float fqrcpar;
  int ncomponent;
  int nsample;
  int nbin;
  int nmode;
  int best2RID;

  float fWeight;
  float rfgWeight;
  float oscwgt;
  TGraph          *byEv_maqe_ccqe_gr;
  TGraph          *byEv_pfo_ccqe_gr;
  TGraph          *byEv_ebo_ccqe_gr;
  TGraph          *byEv_ca5_cc1pi_gr;
  TGraph          *byEv_ca5_ncpiz_gr;
  TGraph          *byEv_ca5_ncpipm_gr;
  TGraph          *byEv_manff_cc1pi_gr;
  TGraph          *byEv_manff_ncpiz_gr;
  TGraph          *byEv_manff_ncpipm_gr;
  TGraph          *byEv_bgscl_cc1pi_gr;
  TGraph          *byEv_bgscl_ncpiz_gr;
  TGraph          *byEv_bgscl_ncpipm_gr;
  TGraph          *byEv_dismpishp_ccoth_gr;
  TGraph          *byEv_sccvec_ccqe_gr;
  TGraph          *byEv_sccvec_ncoth_gr;
  TGraph          *byEv_sccaxl_ccqe_gr;
  TGraph          *byEv_sccaxl_ncoth_gr;
  TGraph          *byEv_rpa_ccqe_gr;
  TBranch          *byEv_maqe_ccqe_br;
  TBranch          *byEv_pfo_ccqe_br;
  TBranch          *byEv_ebo_ccqe_br;
  TBranch          *byEv_ca5_cc1pi_br;
  TBranch          *byEv_ca5_ncpiz_br;
  TBranch          *byEv_ca5_ncpipm_br;
  TBranch          *byEv_manff_cc1pi_br;
  TBranch          *byEv_manff_ncpiz_br;
  TBranch          *byEv_manff_ncpipm_br;
  TBranch          *byEv_bgscl_cc1pi_br;
  TBranch          *byEv_bgscl_ncpiz_br;
  TBranch          *byEv_bgscl_ncpipm_br;
  TBranch          *byEv_dismpishp_ccoth_br;
  TBranch          *byEv_sccvec_ccqe_br;
  TBranch          *byEv_sccvec_ncoth_br;
  TBranch          *byEv_sccaxl_ccqe_br;
  TBranch          *byEv_sccaxl_ncoth_br;
  TBranch          *byEv_rpa_ccqe_br;


  /////////////////////
  //methods
  void setTree(TTree* trin);
  void setTree(TChain* trin);
  void setTree(TTree*, TTree*);
  void setupSplineTree(TTree *h);
  void setParFileName(const char* filename){parFileName=filename;}
  void runPreProcessing();
  void setFVBinHisto();
  void setupNewTree();
  void preProcessIt();
  int passCuts();
  void fillAttributes(t2kfqEvent *fqevent);
  void fillAttributeMap(t2kfqEvent *fqevent);
  int absmode;
  int getComponent();
  int getSample();
  int getBin();
  int getMode();
  int getBest2RFitID();
  float getBANFFWeight();
  float getFixedWeight();
  void processFile(const char* fname,const char* outname="");
  void processFile(const char* f1, const char* f2, const char* outname);
  void processAllFiles(TChain* chain);
  void processAllFiles(TChain*, TChain*);

 private:
  int flgAddMoreVars;
  int flgUseSpikeMask;

};

#endif
