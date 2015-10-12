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



class sift{
  public:
  sift(TChain* chin);
  sift(TTree* trin);
  TTree* tr;
  TTree* trout;
  fqReader* fq;
  visRing*  vis;
  void setupNewTree();
  void siftIt(const char* filename);
  int passCuts();
  int ncomponent;
  int nsample;
  int nbin;
  //new leaf variables
  float towall;
  float wall;
  //selections
  int absmode;
  int getComponent();
  int getSample();
  int getBin();
};

#endif
