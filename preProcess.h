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




using namespace std;


#define NFILEMAX 5000

class preProcess{
  public:
  preProcess(const char* name);
  preProcess(TChain* chin,const char* name="");
  preProcess(TTree* trin,const char* name="");
  void setTree(TTree* trin);
  void setTree(TChain* trin);
  TTree* tr;
  TTree* trout;
  TFile* fout;
  TString nameTag;
  TString fileNames[NFILEMAX];
  int nFiles;
  fqReader* fq;
  visRing*  vis;
  void setupNewTree();
  void preProcessIt();
  int passCuts();
  int ncomponent;
  int nsample;
  int nbin;
  //new leaf variables
  float towall;
  float wall;
  float evtweight;
  float attribute[1000];
  //selections
  void fillAttributes();
  int absmode;
  int getComponent();
  int getSample();
  int getBin();
  float getWeight();
  void processFile(const char* fname,const char* outname="");
  void processAllFiles(TChain* chain);
  void addFile(const char* filename);
  TString getFileRootName();

};

#endif