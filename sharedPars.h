#ifndef __SHAREDPARS_H__
#define __SHAREDPARS_H__

#include <iostream>
#include "TString.h"

#include "keyread.h"

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
  TString preProcessFilesMC;
  TString preProcessFilesData;
  TString preProcessFilesBANFF;
  TString preProcessFilesSpline;
  TString preProcessOutDir;
  int preProcessMCComponents;
  int preProcessFVBinning; 
  int PreProcFCCut;
  double PreProcEVisCut;
  double PreProcWallMinCut;
  double PreProcToWallMinCut;
  int    PreProcNseMax0;
  int    PreProcNseMin;
  double  PreProcInGateCut;
  TString FVBinName0;
  TString FVBinName1;
  TString FVBinName2;
  TString FVBinName3;
  TString FVBinName4;
  TString FVBinName5;
  TString fQAttName0;
  TString fQAttName1;
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
};

#endif
