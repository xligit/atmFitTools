#ifndef HISTCOMPARE_H
#define HISTCOMPARE_H

#include "histoManager.cxx"
#include "TMath.h"
#include "TRandom2.h"
#include "histoTransforms.cxx"
#include "TFitter.h"

using namespace std;


//class to compare histograms and evaluate likelihoods
class histoCompare{
  public:

  //constructors//
  histoCompare(const char* thename);  //standard constructor

  //internal variables
  TString nameTag;  //name associated with this instance
  int nSamp;  //number of samples
  int nBin;  //number of bins
  int nComp;  //number of  components
  int nAtt;  //nummboer of attributes
  //tools for histogram manager management
  //created histo manager from file
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp, int natt);
  histoManager* hManager;
  float Norm;
  float Par[NBINMAX][NCOMPMAX][NATTMAX][2];
  float fixPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  float bestPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  float errPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString parName[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString binName[NBINMAX];
  TString compName[NCOMPMAX];
  TString attName[NCOMPMAX];
  void setBinName(int ibin, const char* name){binName[ibin]=name;}
  void setCompName(int icomp, const char* name){compName[icomp]=name;}
  void setAttName(int iatt, const char* name){attName[iatt]=name;}
  void showFitHisto(int isamp,int ibin,int icomp,int iatt);
  void showFitEffect(int isamp,int ibin,int icomp,int iatt);
  void showFitResult(int isamp,int ibin,int iatt);
  void showFitPars(int ibin,int iatt,int imod);
 // TH1F* hMod[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];

  //tools for adding histogram directly..for debugging
  void addHistogram(TH1F* h,int dataflg);
  TH1F* hData[10];
  TH1F* hMC[10];
  TH1F* hModDebug;
  TH1F* hMod;
  TH1F* hPar;
  TH1F* showSmear(TH1F* h, float smear, float bias);
  void showMod(int imchist);
  TH1F* hTot;
  int nDataHist;
  int nMCHist;
  float parDebug[10][2];
//  float getLDebug(); 
  
  //likelihood evaluateions
  float getSumSq(TH1F* h1, TH1F* h2);
  float getLnL(TH1F* h1, TH1F* h2);
  float getNDiff();
  float getTotSumSq();
  float getTotLnL();
  static void sumSqWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  //for debuggint and play
  float getTotSumSqDebug();
  static void sumSqDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void nDiffDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void minSumSqDebug();
  void minSumSq();
  void sumSqPrefit();
  void LnLFit();
  void LnLPreFit();
  void drawResult(int ihist);

  //staticthis for fits
  static histoCompare* staticthis;
};

#endif
