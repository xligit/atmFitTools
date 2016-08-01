#ifndef SHAREDPARS_C
#define SHAREDPARS_C

#include "TString.h"

#include "keyread.cxx"

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
  void readParsFromFile();

  //shared variables
  TString globalRootName;

  //can get parameters directly by name
  int getParI(const char*);
  double getParD(const char*);
  TString getParS(const char*); 

  //key reading object
  keyread* kr;

  //fit parameters
  int nFVBins;
  int nSamples;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nModes;
  double normFactor;
  TString preProcessFilesMC;
  TString preProcessFilesSpline;
  TString preProcessFilesBANFF;
  TString preProcessFilesData;
  TString preProcessOutDir;
  int preProcessMCComponents;
  int preProcessMCSamples;
  int preProcessFVBinning; 
  int preProcFCCut;
  double preProcEVisCut;
  double preProcWallMinCut;
  double preProcToWallMinCut;
  int    preProcNseMax0;
  int    preProcNseMin;
  int    preProcMaskFlg;
  TString preProcMaskFile;
  TString ntupleType;
  double  preProcInGateCut;
  int     preProcAddMoreVars;
  TString FVBinName0;
  TString FVBinName1;
  TString FVBinName2;
  TString fQAttName0;
  TString fQAttName1;
  TString fQAttName2;
  TString fQAttName3;
  TString fQAttName4;
  TString fQAttName5;
  TString fQAttName6;
  TString fQAttName7;
  TString fQAttName8;
  TString MCComponentName0;
  TString MCComponentName1;
  TString MCComponentName2;
  TString MCComponentName3;
  TString MCComponentName4;
  TString MCComponentName5;
  TString MCComponentName6;
  TString sampleName0;
  TString sampleName1;
  TString sampleName2;
  TString sysParName0;
  TString sysParName1;
  TString sysParName2;
  TString sysParName3;
  TString sysParName4;
  TString sysParName5;
  TString sysParName6;
  TString sysParName7;
  TString sysParName8;
  TString hFactoryOutput;
  TString hFactoryMCFiles;
  TString hFactoryDataFiles;
  TString splineFactoryOutput;
  int MCMCNSteps;
  double MCMCTunePar;
  int useSplinesFlg;
  int fixAllSmearFlg;
  int NMCMCPts;
  int MCMCBurnIn;
  int NMCEvents;
  TString MCMCFile;
  TString sysParType;
  int NDataEvents;
  int flgUseNormPars;
 
};


double sharedPars::getParD(const char* parname){
  return kr->getKeyD(parname);
}

TString sharedPars::getParS(const char* parname){
  return kr->getKeyS(parname);
}

int sharedPars::getParI(const char* parname){
  return kr->getKeyI(parname);
}

void sharedPars::readParsFromFile(){

  //create object to read in keys from file
  kr = new keyread(parFileName.Data());

  //read in contents
  kr->readFile();

  //set parameters to values
  useSplinesFlg = kr->getKeyI("useSplinesFlg");
  MCMCNSteps= kr->getKeyI("MCMCNSteps");
  MCMCTunePar = kr->getKeyD("MCMCTunePar");
  normFactor = kr->getKeyD("normFactor");
  fixAllSmearFlg = kr->getKeyI("fixAllSmearFlg");
  FVBinName0= kr->getKeyS("FVBinName0");
  FVBinName1=kr->getKeyS("FVBinName1");
  FVBinName2=kr->getKeyS("FVBinName2");
  fQAttName0=kr->getKeyS("fQAttName0");
  fQAttName1=kr->getKeyS("fQAttName1");
  fQAttName2=kr->getKeyS("fQAttName2");
  fQAttName3=kr->getKeyS("fQAttName3");
  fQAttName4=kr->getKeyS("fQAttName4");
  fQAttName5=kr->getKeyS("fQAttName5");
  fQAttName6=kr->getKeyS("fQAttName6");
  fQAttName7=kr->getKeyS("fQAttName7");
  fQAttName8=kr->getKeyS("fQAttName8");
  hFactoryDataFiles=kr->getKeyS("hFactoryDataFiles");
  hFactoryMCFiles=kr->getKeyS("hFactoryMCFiles");
  hFactoryOutput=kr->getKeyS("hFactoryOutput");
  MCComponentName0=kr->getKeyS("MCComponentName0");
  MCComponentName1=kr->getKeyS("MCComponentName1");
  MCComponentName2=kr->getKeyS("MCComponentName2");
  MCComponentName3=kr->getKeyS("MCComponentName3");
  MCComponentName4=kr->getKeyS("MCComponentName4");
  MCComponentName5=kr->getKeyS("MCComponentName5");
  MCComponentName6=kr->getKeyS("MCComponentName6");
  sampleName0=kr->getKeyS("sampleName0");
  sampleName1=kr->getKeyS("sampleName1");
  sampleName2=kr->getKeyS("sampleName2");
  sysParName0=kr->getKeyS("sysParName0");
  sysParName1=kr->getKeyS("sysParName1");
  sysParName2=kr->getKeyS("sysParName2");
  sysParName3=kr->getKeyS("sysParName3");
  sysParName4=kr->getKeyS("sysParName4");
  sysParName5=kr->getKeyS("sysParName5");
  sysParName6=kr->getKeyS("sysParName6");
  sysParName7=kr->getKeyS("sysParName7");
  sysParName8=kr->getKeyS("sysParName8");
  nFVBins     = kr->getKeyI("nFVBins");
  nSamples    = kr->getKeyI("nSamples");
  nComponents = kr->getKeyI("nComponents");
  nAttributes  = kr->getKeyI("nAttributes"); 
  nSysPars    = kr->getKeyI("nSysPars");
#ifdef NMODE
  nModes = NMODE;
#endif
  preProcessFilesMC = kr->getKeyS("preProcessFilesMC"); 
  preProcessOutDir = kr->getKeyS("preProcessOutDir"); 
  preProcessFilesData = kr->getKeyS("preProcessFilesData"); 
  preProcessFilesSpline = kr->getKeyS("preProcessFilesSpline");
  preProcessFilesBANFF = kr->getKeyS("preProcessFilesBANFF");
  preProcessMCComponents = kr->getKeyI("preProcessMCComponents");
  preProcessFVBinning = kr->getKeyI("preProcessFVBinning");
  preProcessMCSamples = kr->getKeyI("preProcessMCSamples");
  preProcAddMoreVars = kr->getKeyI("preProcAddMoreVars");
  preProcMaskFile = kr->getKeyS("preProcMaskFile");
  preProcMaskFlg = kr->getKeyI("preProcMaskFlg");
  ntupleType=kr->getKeyS("ntupleType");
  globalRootName = kr->getKeyS("globalRootName");
  splineFactoryOutput = kr->getKeyS("splineFactoryOutput");
  sysParType = kr->getKeyS("sysParType");
  NMCMCPts = kr->getKeyI("NMCMCPts");
  MCMCBurnIn=kr->getKeyI("MCMCBurnIn");
  NMCEvents=kr->getKeyI("NMCEvents");
  MCMCFile=kr->getKeyS("MCMCFile");
  preProcFCCut=kr->getKeyI("preProcFCCut");;
  preProcEVisCut=kr->getKeyD("preProcEVisCut");
  preProcWallMinCut=kr->getKeyD("preProcWallMinCut");
  preProcToWallMinCut=kr->getKeyD("preProcToWallMinCut");
  preProcNseMax0=kr->getKeyI("preProcNseMax");
  preProcNseMin=kr->getKeyI("preProcNseMin");
  preProcInGateCut=kr->getKeyD("preProcInGateCut");
  NDataEvents = kr->getKeyI("NDataEvents");
  flgUseNormPars = kr->getKeyI("flgUseNormPars");
}

sharedPars::sharedPars(const char* parfilename){
  parFileName = parfilename;
//  readParsFromFile();
}


#endif
