#ifndef HISTCOMPARE_H
#define HISTCOMPARE_H

#include "histoManager.C"
#include "TMath.h"
#include "TRandom2.h"
#include "histoTransforms.C"
#include "TFitter.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "likelihood.C"

using namespace std;


//class to compare histograms and evaluate likelihoods
class histoCompare{
  public:

  //constructors//
  histoCompare(const char* parfile);  //construct from parameter file
  histoCompare();  //standard constructor

  //internal variables
  TString nameTag;  //name associated with this instance
  int nSamp;  //number of samples
  int nBin;  //number of bins
  int nComp;  //number of  components
  int nAtt;  //nummboer of attributes
  double tunePar; //tuning parameter for MCMC
  //tools for histogram manager management
  //created histo manager from file
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp, int natt);
  histoManager* hManager;
  //atmospheric pars
  atmFitPars* thePars;
  sharedPars* runPars;
  int MCMCNSteps;
  double Norm;
  double Par[NBINMAX][NCOMPMAX][NATTMAX][2];
  double sysPar[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  
  double fixPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  double bestPar[NBINMAX][NCOMPMAX][NATTMAX][2];
 // double errParLo[NBINMAX][NCOMPMAX][NATTMAX][2];
 // double errParHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  double errParLo[1000];
  double errParHi[1000];

  TString parName[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString binName[NBINMAX];
  TString compName[NCOMPMAX];
  TString attName[NCOMPMAX];
  void setBinName(int ibin, const char* name){binName[ibin]=name;}
  void setCompName(int icomp, const char* name){compName[icomp]=name;}
  void setAttName(int iatt, const char* name){attName[iatt]=name;}
  void setupPars(int nsyspars=0); //sets up all parameters
  //post-fit toolts
  void profileL(int ibin, int icomp, int iatt, int imod, double range, int npts=100);
  void profileL(int ipar,double range, int npts=100,int sameflg=0);
  void showFitHisto(int isamp,int ibin,int icomp,int iatt);
  void showFitEffect(int isamp,int ibin,int icomp,int iatt);
  void showFitResult(int isamp,int ibin,int iatt);
  void showFitPars(int ibin,int iatt,int imod);
  void showModHiso(int isamp,int ibin, int icomp, int iatt, double smear, double bias);
  void runMCMC(int nsteps);
 // double getErrLo(int ibin,int icomp,int iatt,int imod);
  double getErrLo(int isyst);
  double getErrHi(int isyst);
  TH1D* getModifiedHisto(int ibin, int icomp, int iatt){return hManager->getSumHistogramMod(ibin,icomp,iatt);}
 // TH1D* hMod[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];

  //initialize all necessary components
  void initialize(histoManager* hm, atmFitPars* apars);

  void setupSplines(const char* fname,int npar){hManager->readSplinesFromFile(fname,npar);}
  //tools for adding histogram directly..for debugging
  void addHistogram(TH1D* h,int dataflg);
  int  rebinFactor;
  void setRebinFactor(int ifact){rebinFactor=ifact;}
  TH1D* hData[10];
  TH1D* hMC[10];
  TH1D* hModDebug;
  TH1D* hMod;
  TH1D* hTmp; //temporary histogram container
  TH1D* hPar; //for fit parameters
  TGraphAsymmErrors* gPar;
  TH1D* hParErrLo;
  TH1D* hParErrHi;
  TH1D* hProf[2];
  TH1D* showSmear(TH1D* h, double smear, double bias);
  void showMod(int imchist);
  TH1D* hTot;
  TCanvas* cc;
  int nDataHist;
  int nMCHist;
  int useLnLType;
  double parDebug[10][2];
//  double getLDebug(); 
 
  double cScale; //correction scale for likelihood
  void printFitResults(const char* directory);

  //////////////////////////////////////////////////////////
  //flags
  int flgFixAllSmearPars;
 
  //likelihood evaluateions
  double getSumSq(TH1D* h1, TH1D* h2);
  double getLnL(TH1D* h1, TH1D* h2);
  double getNDiff();
  double getTotSumSq();
  double getTotLnL();
  static void sumSqWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void  getTotLnL1D(double& result,int npar, double par[]);
  //for debuggint and play
  double getTotSumSqDebug();
  static void sumSqDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void nDiffDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void minSumSqDebug();
  void minSumSq();
  void sumSqPrefit();
  void LnLFit();
  void LnLPreFit();
  void singleParFit(int ipar);
  void sysParFit();
  void drawResult(int ihist);
  void timetest(int ntry);
  void saveFitPars(const char* filename); //< write parameters and outputs to a file
  void readFitPars(const char* filename); //< read parameters from a file
  void tuneMCMC(int ncyles=1,int nsteps=150,double goal=0.25);
  void tuneMCMCOld(int ncyles=1,int nsteps=150,double goal=0.25);
  void calcRoughParErr();
  //staticthis for fits
  static histoCompare* staticthis;
};

#endif
