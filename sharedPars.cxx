#include "sharedPars.h"

double sharedPars::getParD(const char* parname){
  return kr->getKeyD(parname);
}

TString sharedPars::getParS(const char* parname){
  return kr->getKeyS(parname);
}

int sharedPars::getParI(const char* parname){
  return kr->getKeyI(parname);
}

void sharedPars::readParsFromFile(const char* filename){
  //if file name is given as argument, use it
//  if (filename){
//    cout<<"sharedPars: setting input file name to : "<<filename<<endl;
//    parFileName = filename;
//  }
  
  //create object to read in keys from file
  kr = new keyread(parFileName.Data());

  //read in contents
  kr->readFile();

  //set parameters to values
  useSplinesFlg = kr->getKeyI("useSplinesFlg");
  MCMCNSteps= kr->getKeyI("MCMCNSteps");
  MCMCTunePar = kr->getKeyD("MCMCTunePar");
  fixAllSmearFlg = kr->getKeyI("fixAllSmearFlg");
  FVBinName0= kr->getKeyS("FVBinName0");
  FVBinName1=kr->getKeyS("FVBinName1");
  FVBinName2=kr->getKeyS("FVBinName2");
  FVBinName3=kr->getKeyS("FVBinName3");
  FVBinName4=kr->getKeyS("FVBinName4");
  FVBinName5=kr->getKeyS("FVBinName5");
  fQAttName0=kr->getKeyS("fQAttName0");
  fQAttName1=kr->getKeyS("fQAttName1");
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
  preProcessFilesMC = kr->getKeyS("preProcessFilesMC"); 
  preProcessOutDir = kr->getKeyS("preProcessOutDir"); 
  preProcessFilesData = kr->getKeyS("preProcessFilesData"); 
  preProcessFilesBANFF = kr->getKeyS("preProcessFilesBANFF");
  preProcessFilesSpline = kr->getKeyS("preProcessFilesSpline");
  preProcessMCComponents = kr->getKeyI("preProcessMCComponents");
  preProcessFVBinning = kr->getKeyI("preProcessFVBinning");
  globalRootName = kr->getKeyS("globalRootName");
  splineFactoryOutput = kr->getKeyS("splineFactoryOutput");
  sysParType = kr->getKeyS("sysParType");
  NMCMCPts = kr->getKeyI("NMCMCPts");
  MCMCBurnIn=kr->getKeyI("MCMCBurnIn");
  NMCEvents=kr->getKeyI("NMCEvents");
  MCMCFile=kr->getKeyS("MCMCFile");
  PreProcFCCut=kr->getKeyI("PreProcFCCut");;
  PreProcEVisCut=kr->getKeyD("PreProcEVisCut");
  PreProcWallMinCut=kr->getKeyD("PreProcWallMinCut");
  PreProcToWallMinCut=kr->getKeyD("PreProcToWallMinCut");
  PreProcNseMax0=kr->getKeyI("PreProcNseMax");
  PreProcNseMin=kr->getKeyI("PreProcNseMin");
  PreProcInGateCut=kr->getKeyD("PreProcInGateCut");
}

sharedPars::sharedPars(const char* parfilename){
  parFileName = parfilename;
//  readParsFromFile();
}
