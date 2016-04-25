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
#include "TGraph.h"
#include <map>
#include <string>
#include "TH2FV.C"
using namespace std;


#define NFILEMAX 5000


//////////////////////////////////////////////////////////////
//Class to take care of preprocessing of all data and MC files
//Usage:
//	1) create instance using preProcess()
//	2) specify parameter file using setParFileName(<name>)
//      3) run using runPreProcessing()
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
  int MCSamples; //< flag for sample definitions
  TH1D* hWeight;
  TGraph* gWeight;
  int useWeights;
  map<string,double> attributeMap;
  TString attributeList[50];
  int nAttributes;
  TH2FV* hFVBins;

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
  int best2RID;

  
  /////////////////////
  //methods
  void setTree(TTree* trin);
  void setTree(TChain* trin);
  void setParFileName(const char* filename){parFileName=filename;}
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
  

};

#endif
