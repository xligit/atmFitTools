#ifndef HISTOFACTORY_H
#define HISTOFACTORY_H

#include "histoManager.C"
//#include "shared.h"

using namespace std;

//class for building histograms from sorted data and mc trees
class histoFactory{
  public:

  histoFactory(int nsampl,int nbins,int ncomp,const char* name=""); //creates blank histogram factory
  histoManager* hManager; 
  TTree* dataTree; //tree containing data events
  TTree* mcTree;  //tree conttaining MC events
  fQreader* fqData;  //reads data tree
  fQreader* fqMC; //reads MC tree
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TH1F*    htmp; //temporary histogram pointer for drawing and comparison
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  float att[NATTMAX]; //array of all attribute values
  TString attNames[NATTMAX];  //array of attribute names
  TString attType[NATTMAX];  //array of attribute type codes
  void init();  //initialize after attributes have been set (sets branch addresses, creates histograms)
  void addAttribute(int iatt);  //add an attribute (fiTQun variable) to list of histograms to be made
  TH1F* getHistogram(int iatt,const char* thename); //returns pointer to MC histogram
  TH1F* getHistogramData(int iatt,const char* thename); //returns pointer to Data histogram
  void fillAttributesData(); //fills all data histograms
  void fillAttributesMC();  //fills all MC histograms
  void fillHistos(); //fills all histograms
  //setters
  void setDataTree(TTree* tr);
  void setDataTree(TChain* ch);
  void setMCTree(TTree* tr);
  void setMCTree(TChain* ch);
  void saveToFile();
  TString outputFileName;
  TString getOutputFileName(){return outputFileName;}
};


#endif
