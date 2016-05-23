#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__

//#include "/nfs/hepusers/users/amissert/stdROOTinc.h"
#include <iostream>
#include "visRing.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph.h"
#include "TString.h"
#include "FVCalculators.cxx"
#include "sharedPars.h"
#include "shared.h"

using namespace std;

#define NFILEMAX 5000

//////////////////////////////////////////////////////////////
//Class to take care of preprocessing of all data and MC files
//Usage:
//  1) Specify parameters in sharedpars.dat file.  These parameters 
//     specify the files to be run on, where the output should go,
//     and which FV binning and MC component definitions to use
//  2) Call runPreprocess()

class preProcess{
  public:

  /////////////////////
  //constructors
  preProcess();
  preProcess(TChain* chin,const char* name="");
  preProcess(TTree* trin,const char* name="");
  preProcess(TChain *mc, TChain *spline, const std::string name);
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
  fqReader* fq;
  visRing*  vis;
  int MCComponents;
  int FVBinning;
  int nFiles;
  int nFilesSpline;
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

  //////////////////////
  //new branch variables
  float towall;
  float wall;
  float evtweight; // banff weight
  float rfgweight;
  float attribute[1000];
  int ncomponent;
  int nsample;
  int nbin;
  int nmode;
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
  void setupNewTree();
  void preProcessIt();
  int passCuts();
  void fillAttributes();
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

};

#endif
