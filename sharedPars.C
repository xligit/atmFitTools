#include "TString.h"

#include "keyread.C"

#include <iostream>

using namespace std;

//class to hold and read in shared parameters for fits
class sharedPars{
  public:

  //constructor
  sharedPars(const char* parfilename);

  //name of parameter file
  TString parFileName;
  
  //read in values from parameter file
  void readParsFromFile(const char* filename="");

  //shared variables
  TString globalRootName;

  //fit parameters
  int nFVBins;
  int nSamples;
  int nComponents;
  int nAttributes;
  int nSysPars;
  
};

void sharedPars::readParsFromFile(const char* filename){
  //if file name is given as argument, use it
  if (filename!=""){
    cout<<"sharedPars: setting input file name to : "<<filename<<endl;
    parFileName = filename;
  }
  
  //create object to read in keys from file
  keyread* kr = new keyread(parFileName.Data());

  //read in contents
  kr->readFile();

  //set parameters to values
  nFVBins     = kr->getKeyI("nFVBins");
  nSamples    = kr->getKeyI("nSamples");
  nComponents = kr->getKeyI("nComponents");
  nAttributes  = kr->getKeyI("nAttributes"); 
  nSysPars    = kr->getKeyI("nSysPars");
  
  globalRootName = kr->getKeyS("globalRootName");
}

sharedPars::sharedPars(const char* parfilename){
  parFileName = parfilename;
  readParsFromFile();
}

