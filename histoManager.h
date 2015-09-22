#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H

#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "fQreader.C"
#include <iostream>
#include "THStack.h"
#include "TLegend.h"

#define NSAMPMAX 5
#define NCOMPMAX 20
#define NATTMAX 20
#define NBINMAX 10

class histoManager{
  public:
  histoManager(int nsampl,int nbins,int ncomp,const char* name=""); //constructor
  histoManager(const char* rootname,int nsamp,int nbin,int ncomp,int natt);
  TTree* dataTree;
  TTree* mcTree; 
  fQreader* fqData;
  fQreader* fqMC;
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TFile*   fin; //input file of histos to be read
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  TH1F* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TH1F* hData[NSAMPMAX][NBINMAX][NATTMAX];  //array of all Data histograms
  TLegend* Leg;
  float att[NATTMAX]; //array of all attribute values
  TString attNames[NATTMAX];  //array of attribute names
  TString attType[NATTMAX];  //array of attribute type codes
  void init();  //initialize after attributes have been set
  void addAttribute(int iatt);  //add an attribute (fiTQun variable) to list
  
  //histogram maker
  TH1F* getHistogram(int iatt,const char* thename); 
  TH1F* getHistogramData(int iatt,const char* thename); 
  //fill attribute variables
  void fillAttributesData();
  void fillAttributesMC();
  void fillHistos();

  //setters
  void setDataTree(TTree* tr);
  void setDataTree(TChain* ch);
  void setMCTree(TTree* tr);
  void setMCTree(TChain* ch);

  //plotting
  void showMCBreakdown(int isample,int ibin,int iatt);
  THStack* showMCBreakdownStack(int isample,int ibin,int iatt);
  TH1F* calcMCSum(int isample, int ibin, int iatt);
 
  //file management
  void saveToFile();
  void readFromFile(const char* rootename,int nsamp,int nbin,int ncomp,int natt);
};

#endif
