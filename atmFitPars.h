#ifndef __ATMFITPARS_H__
#define __ATMFITPARS_H__

#include "shared.h"
#include "TRandom2.h"
#include "sharedPars.h"
#include "covXsec.h"
#include "covBANFF.h"
#include "covBase.h"

using namespace std;

//class containing all atmospheric fit pars
class atmFitPars{
  public:
  
  ///////////////////////////////////////////////////////////////////
  //constructors
  atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst=0);
  atmFitPars(int isamp, int ibin, int icomp, int iatt, const char* systype); 
  atmFitPars(const char* parfile, covBase *covm = 0); //constructs from parameter file

  ///////////////////////////////////////////////////////////////
  //numbers of various parametrs
  int nSamples; // # of sub-events
  int nBins; // fiducial volume bin
  int nComponents; // 1Re, 1Rmu, MR1e, MR1mu, 2Rpi0, other
  int nAttributes; // pid, ring-counting
  int nSysPars; 
  int nTotPars;
  int nModes;
  sharedPars* runpars;
  covBase *cov;
  TRandom3 *rnd;
  ///////////////////////////////////////////////////////////////////
  //parameter values
  double histoPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncLo[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  double sysPar[NSYSPARMAX];
  double sysParNom[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  double sysParUp[NSYSPARMAX];
  double sysParLow[NSYSPARMAX];
  std::string sysParName[NSYSPARMAX];
  double pars[4000]; // current values
  double parsProp[4000];
  double parUnc[4000];
  int   fixPar[4000]; //< array of fix flags for parameters
  double bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2];
  double norm;  
  float fScale;
  //////////////////////////////////////////////////////////////
  //methods
  void proposeStep();
  void acceptStep();
  void setStepSize(float f) {fScale = f; cov->setStepScales(f);}
  void setSeed(int i) {rnd->SetSeed(i);}
  void setNorm(double x) {norm=x;}
  void initPars(const char* systype=""); //< sets parameters to initial values
  int getParIndex(int ibin, int icomp, int iatt, int imod){return parIndex[ibin][icomp][iatt][imod];}
  double getParameter(int ipar){return pars[ipar];}
  double getPropParameter(int ipar) { return parsProp[ipar]; }
  double getSysParameter(int ipar) { return sysPar[ipar]; }
  double* getParameters() {return pars;}
  void setParameter(int ipar, double value);
  void setSysParameter(int ipar, double value);
  void setParameter(int ibin, int icomp, int iatt, int imod, double value); 
  void setSysParUnc(int isys,double value) {sysParUnc[isys]=value;}
  void fixParameter(int ipar);
  void fixParameter(int ibin,int icomp,int iatt, int imod);
  void fixAllSmearPars(int isfixed=1);
  void setRandSysPar(); //sets systematic parameters to random values
  int  checkFixFlg(int ibin,int icomp,int iatt, int imod);
  void resetDefaults();
  void printParValues();
  void setCov(covBase *covariance);
  //void setSysParSig(double sigma = 0.);

  int binOfPar[4000];
  int compOfPar[4000];
  int attOfPar[4000];
  int typeOfPar[4000];
  std::string sysType;

  //saving and reading pars
  void savePars(const char* filename);
  void readPars(const char* filename);
  void printPars();
};

#endif
