#ifndef SIFT_H
#define SIFT_H

//#include "/nfs/hepusers/users/amissert/stdROOTinc.h"
#include <iostream>
#include "visRing.C"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "FVCalculators.C"
#include "sharedPars.C"



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


  /////////////////////
  //internal variables
  TChain* chmc;
  TChain* chdat;
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
  float evtweight;
  float attribute[1000];
  int ncomponent;
  int nsample;
  int nbin;

  
  /////////////////////
  //methods
  void setTree(TTree* trin);
  void setTree(TChain* trin);
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
  int getBest2RFitID();
  float getWeight();
  void processFile(const char* fname,const char* outname="");
  void processAllFiles(TChain* chain);

};

#endif
