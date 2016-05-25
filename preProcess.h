#ifndef PREPROCESS_H
#define PREPROCESS_H

//#include "/nfs/hepusers/users/amissert/stdROOTinc.h"
#include <iostream>
#include "visRing.C"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "FVCalculators.C"
#include "sharedPars.C"
#include "TGraph.h"
#include <map>
#include <string>
#include "TH2FV.C"
#include "masktools.C"
#include "TObjArray.h"

#define NFILEMAX 5000

using namespace std;



//////////////////////////////////////////////////////////////
// Class to take care of preprocessing of all data and MC files
// Usage:
//	1) create instance using preProcess()
//	2) specify parameter file using setParFileName(<name>)
//  3) run using runPreProcessing()

class preProcess{
  public:

  /////////////////////
  // constructors
  preProcess();


  /////////////////////
  // internal variables
  TChain* chmc; //< points to input MC
  TChain* chdat; //< points to input data
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
  int MCSamples; //< flag for sample definitions
  TH1D* hWeight;
  TGraph* gWeight;
  int useWeights;
  map<string,double> attributeMap;
  TString attributeList[50];
  int nAttributes;
  TH2FV* hFVBins;
  TH1D* hmask;

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
  float towallv[50];
  float minconev[50];
  float perimv[50];
  float wallv2;
  float fq1rwall[10][7];
  float fq1rtowall[10][7];
  float fq1rmincone[10][7];
  float fq1rperim[10][7];
  float evtweight;
  float attribute[1000];
  int ncomponent;
  int nsample;
  int nbin;
  int best2RID;

  
  /////////////////////
  //methods
  void setTree(TTree* trin);
  void setTree(TChain* trin);
  void setParFileName(const char* filename);
  void runPreProcessing();
  void setFVBinHisto();
  void setupNewTree();
  int preProcessIt(); //< process a file and return the number of entries in new tree
  int passCuts();
  void fillAttributes(fqReader* fqevent);
  void fillAttributeMap(fqReader* fqevent);
  int absmode;
  int getComponent();
  int getSample();
  int getBin();
  int getBest2RFitID(); //< returns best 2R fit ID for current event
  float getWeight();
  void processFile(const char* fname,const char* outname="");
  void processAllFiles(TChain* chain);
  //sets a histogram to calculate weights for events (use to correct cosmic
  //muon momenum distribution)
  void setWeightHistogram(const char* file, const char * name);
 
  private:

  int flgAddMoreVars;
  int flgUseSpikeMask;
};

#endif

#ifndef PREPROCESS_C
#include "preProcess.C"
#endif
