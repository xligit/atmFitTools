#ifndef HISTCOMPARE_H
#define HISTCOMPARE_H

#include "histoManager.cxx"
#include "TMath.h"
#include "TRandom2.h"
#include "histoTransforms.cxx"
#include "TFitter.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "likelihood.cxx"
#include "TGraph2D.h"
#include "TRandom3.h"
#include "shared.h"

using namespace std;

//class to compare histograms and evaluate likelihoods
class histoCompare{
  public:

  //constructors//
  histoCompare(const char* parfile, bool separateneutmode = false);  //construct from parameter file
  histoCompare();  //standard constructor

  //internal variables
  bool separateNeutMode;
  TString nameTag;  //name associated with this instance
  int nSamp;  //number of samples
  int nBin;  //number of bins
  int nComp;  //number of  components
  int nAtt;  //nummboer of attributes
  int nMode;
  double tunePar; //tuning parameter for MCMC
  //tools for histogram manager management
  //created histo manager from file
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp, int natt);
  // For T2K parameterization //
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp, int imode, int natt);
  //
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
  TString diffChainFileName;

  TString parName[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString binName[NBINMAX];
  TString compName[NCOMPMAX];
  TString attName[NCOMPMAX];
  void setBinName(int ibin, const char* name){binName[ibin]=name;}
  void setCompName(int icomp, const char* name){compName[icomp]=name;}
  void setAttName(int iatt, const char* name){attName[iatt]=name;}
  // For T2K parameterization
  void setupPars(int nsyspars=0); //sets up all parameters  
  void setupPars(atmFitPars *a);

 
  //
  //post-fit tools and plotting
  void profileL(int ibin, int icomp, int iatt, int imod, double range, int npts=100);
  void profileL(int ipar,double range, int npts=100,int sameflg=0);
  void showFitHisto(int isamp,int ibin,int icomp,int iatt);
  void showFitEffect(int isamp,int ibin,int icomp,int iatt);
  void showFitResult(int isamp,int ibin,int iatt);
  void showFitDiff(int isamp, int ibin, int iatt);
  void showFitPars(int ibin,int iatt,int imod);
  void setPar(int ipar, double value){thePars->setParameter(ipar,value);}
  void showModHiso(int isamp,int ibin, int icomp, int iatt, double smear, double bias);
  TH2D* show2DLnL(int parx, double xmin, double xmax, int pary, double ymin, double ymax, int npts=100);
  TGraph2D* show2DLnLG(int parx, double xmin, double xmax, int pary, double ymin, double ymax, int npts=100);
  void runMCMC(int nsteps);
  void runDEMCMC(int nstep);
  void runDiffMCMC(int nsteps); //< fill with differential steps
 // double getErrLo(int ibin,int icomp,int iatt,int imod);
  double getErrLo(int isyst);
  double getErrHi(int isyst);
  TH1D* getModifiedHisto(int ibin, int icomp, int iatt){return hManager->getSumHistogramMod(ibin,icomp,iatt);}
 // TH1D* hMod[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];

  //initialize all necessary components
  void initialize(histoManager* hm, atmFitPars* apars);

  void setupSplines(const char* fname){hManager->readSplinesFromFile(fname);}
  //tools for adding histogram directly..for debugging
  void printParName(int ipar){cout<<thePars->parName[ipar].Data()<<endl;}
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
  int flgUsePriorsInFit; 

  /////////////////////////////////////////////////////////
  // make fake data by setting modified histograms as data
  void makeFakeData();

  //likelihood evaluateions;
  double getSumSq(TH1D* h1, TH1D* h2);
  double getLnL(TH1D* h1, TH1D* h2);
  double getTotSumSq();
  double getTotLnL();
  void  getTotLnL1D(double& result,int npar, double par[]);
  void LnLFit();
  void LnLPreFit();
  static void lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void singleParFit(int ipar);
  void sysParFit();
  void timetest(int ntry);
  void saveFitPars(const char* filename); //< write parameters and outputs to a file
  void readFitPars(const char* filename); //< read parameters from a file
  void tuneMCMC(int ncyles=1,int nsteps=150,double goal=0.25);
  void tuneDEMCMC(int ncyles=1,int nsteps=150,double goal=0.25);
  void tuneMCMCOld(int ncyles=1,int nsteps=150,double goal=0.25);
  void calcRoughParErr();

  //staticthis for fits
  static histoCompare* staticthis;
};

#endif
