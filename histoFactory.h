#ifndef HISTOFACTORY_H
#define HISTOFACTORY_H

#include "shared.h"
#include "sharedPars.cxx"
#include "fqProcessedEvent.h"
#include "histoManager.cxx"

using namespace std;

/////////////////////////////////////////////////////////////////////
//class for building histograms from processed data and mc trees
class histoFactory{
  public:
  
  ////////////////////////////////////////////////////////////////////
  // Constructors
  //obsolete, do not use
  histoFactory(int nsampl,int nbins,int ncomp,const char* name=""); //creates blank histogram factory
  // use ths one 
  histoFactory(const char* parfile);//< initialize using parameters in par file

  ///////////////////////////////////////////////////////////////
  //Internal variables
  histoManager* hManager; //manages pointers to hisograms 
  TTree* dataTree; //tree containing data events
  TTree* mcTree;  //tree conttaining MC events
  fqProcessedEvent* fqData;  //reads data tree
  fqProcessedEvent* fqMC; //reads MC tree
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TH1D*    htmp; //temporary histogram pointer for drawing and comparison
  TH1D*    hnorm;
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  int nDataEvents; //number of data events
  int nMCEvents;  //number of MC events
  double normMC; //normalization factor for MC
  double att[NATTMAX]; //array of all attribute values
  TString parFileName; //name of parameter file at initialization
  sharedPars* runpars; //runtime paramters from parameter file

  /////////////////////////////////////////////////////////////
  //Methods
  void init();  //initialize after attributes have been set (sets branch addresses, creates histograms)
  TH1D* getHistogram(int iatt,const char* thename); //returns pointer to MC histogram
  TH1D* getHistogramData(int iatt,const char* thename); //returns pointer to Data histogram
  void fillHistos(); //fills all histograms
  void setDataTree(TTree* tr);
  void setDataTree(TChain* ch);
  void setMCTree(TTree* tr);
  void setMCTree(TChain* ch);
  void saveToFile();
  TString outputFileName;
  TString getOutputFileName(){return outputFileName;}
  void setOutputFileName(const char *name){outputFileName=name;}
  void normalizeHistos(double scale=-1.);
  void runHistoFactory(); //< this will run the factory and fill all histograms
};

#endif

#ifndef HISTOFACTORY_C
#include "histoFactory.cxx"
#endif

