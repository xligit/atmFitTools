#ifndef __T2KSPLINEFACTORY_H__
#define __T2KSPLINEFACTORY_H__

#include "hSplines.cxx"
#include "histoManager.h"
#include "sharedPars.cxx"
#include "atmFitPars.h"
#include "shared.h"
#include "t2kfqReader.h"

//#define NSYSTMAX 500
//#define NHBINSMAX 300
//#define NPTSMAX 21 // used to be 5    
//#define NMODE 9

//class for creating splines
class t2kSplineFactory{
  public:
  //constructor
  t2kSplineFactory(int nsamp, int nbin, int ncomp, int natt, int nsyst, const char* name, int nmode = 0, bool separateneutmode = false);
  //t2kSplineFactory(int nsamp, int nbin, int ncomp, int nmode, int natt, int nsyst, const char *name);
  t2kSplineFactory(const char* parfile, bool separateneutmode = false);//< initialize using parameters in par file
  //t2kSplineFactory();
  
  //internal vars
  TString parFileName; //< name of parameter file
  TString nameTag; //< set in constructor. this is the prefix for the output file
  TChain* mcTree; 
  t2kfqReader* mcEvt;
  histoManager* hManager; //manages all default histograms
  atmFitPars *fitPars;
  TH1D* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; //array for modified histograms for spline creation
  TH1D *hMCMode0[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX];
  TH1D *hMCMode1[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode2[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode3[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode4[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode5[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode6[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode7[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode8[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode9[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX];
  void resetModHistos();
  TFile *fout; //output file
  TString foutName;
  TH1D* htmp; //temporary histogram pointer
  int nSamp;
  int nBin;
  int nComp;
  int nAtt;
  int nMode;
  int nSyst;
  bool separateNeutMode;
  TString sysParType; //< code denoting the type of parameterization used, see setupSysPar
  double sysPar[NSYSPARMAX]; //systematic parameter values
  double sysUnc[NSYSPARMAX];  //systematic parameter uncertainties
  std::string sysName[NSYSPARMAX];
  int nP[NSYSPARMAX];
  double attribute[NATTMAX];
  double eventWeight;
  sharedPars* runpars; //< runtime parameters

  //for output tree
  TTree* splineTree;
  int nbin;
  int ncomponent;
  int nattribute;
  int nsample;
  int nmode;
  int nsystpar;
  int npoints;
  int nhistobins;
  double systParValues[NPTSMAX];
  double binWeight[NPTSMAX][NHBINSMAX];

  const static double sigvals[13];
  const static double maqevals[21];
  const static double binarys[2];
  //methods
  //this needs to be modified for each systematic paramater to add
  double getEvtWeight(int ipar); //returns event weight after applying syst. par. 
  double getEvtWeight(t2kfqReader* mcEvt,int ipar,double value); //
  void setOutputFileName(const char* name){foutName=name;}
  TString getOutputFileName(){return foutName;}
  //
  void fillHistograms(int ipt,int isyst); //fills all histograms given weight
  void  makeManagerFromFile(const char* fname); //reads in histograms from histoFactory
  void fillLeaves(int nsamp,int nbin,int ncomp,int natt,int isyst, int imode = -1); //fills leaves of output tree
  void setMCTree(TChain* tr);
  //build the splines
  void buildTheSplines();

  //debugging
  void debugtest();

  //do everything
  void runSplineFactory();

//  private:
  void fillAttributes();
  int getBest2RFitID();
  void setupHistos();
  void setupSystPars(); //sets up systematic parameters
  void incrementSystPars(double nsig);
  void incrementSystPars(double nsig, int i);
  void setAtmFitPars(atmFitPars *a);
};

#endif
