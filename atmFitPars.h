#ifndef ATMFITPARS_H
#define ATMFITPARS_H

#include "shared.h"
#include "TRandom2.h"
#include "sharedPars.cxx"

#ifdef T2K
#include "covXsec.h"
#include "covBANFF.h"
#include "covBase.h"
#endif

using namespace std;

//class containing all atmospheric fit pars
class atmFitPars{
  public:
  
  ///////////////////////////////////////////////////////////////////
  //constructors
  atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst=0);
  atmFitPars(int isamp, int ibin, int icomp, int iatt, const char* systype); 
  // use this constructor
  atmFitPars(const char* parfile); //constructs from parameter file
#ifdef T2K
  atmFitPars(const std::string parfile, covBase *covm = 0);
#endif

  ///////////////////////////////////////////////////////////////
  //numbers of various parametrs
  int nSamples; // # of sub-events
  int nBins; // fiducial volume (and evis) bin 
  int nComponents; // 1Re, 1Rmu, MR1e, MR1mu, 1Rpi0, other 
  int nAttributes; // pid, ring-counting, etc
  int nSysPars;
  int nModes;
  int nTotPars;
  int flgUseNormPars;
  sharedPars* runpars;
#ifdef T2K
  covBase *cov;
  TRandom3 *rnd;
#endif

  ///////////////////////////////////////////////////////////////////
  //parameter values
  double histoNorm[NSAMPMAX][NBINMAX];
  double histoPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncLo[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  double sysPar[NSYSPARMAX];
  double sysParDefault[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  double pars[4000];
  double parDefaultValue[4000];
  double parUnc[4000];
  int   fixPar[4000]; //< array of fix flags for parameters
  double bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2]; //< stores 1D array position for bias/smear pars
  int   sysParIndex[NSYSPARMAX]; //< stores 1D array position for systematic pars
  int   normParIndex[NSAMPMAX][NBINMAX]; //< stores 1D array position for normalization pars
  double norm;  
  // for T2K parameterization///////////////////////////////////
  float fScale;
  double sysParNom[NSYSPARMAX]; // NEUT "nominal" xsec parameter value
  double sysParUp[NSYSPARMAX];
  double sysParLow[NSYSPARMAX];
  std::string sysParName[NSYSPARMAX];
  double parsProp[4000];
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  //methods
  void setNorm(double x){norm=x;}
  void initPars(const char* systype=""); //< sets parameters to initial values
  int getParIndex(int ibin, int icomp, int iatt, int imod){return parIndex[ibin][icomp][iatt][imod];}
  double getParameter(int ipar){return pars[ipar];}
  double getHistoParameter(int ibin, int icomp, int iatt, int imod);
  double getSysParameter(int isys);
  double getNormParameter(int isamp, int ibin);
  void setParameter(int ipar, double value);
  void setSysParameter(int ipar, double value);
  void setParameter(int ibin, int icomp, int iatt, int imod, double value); 
  void setSysParUnc(int isys,double value){sysParUnc[isys]=value;}
  void fixParameter(int ipar);
  void fixParameter(int ibin,int icomp,int iatt, int imod);
  void fixAllSmearPars(int isfixed=1);
  void setRandSysPar(); //sets systematic parameters to random values
  int  checkFixFlg(int ibin,int icomp,int iatt, int imod);
  void resetDefaults();
  void printParValues();
  int binOfPar[4000];
  int compOfPar[4000];
  int attOfPar[4000];
  int typeOfPar[4000];
#ifdef T2K
  void proposeStep();
  void acceptStep();
  void setStepSize(float f) {fScale = f; cov->setStepScales(f);}
  void setSeed(int i) {rnd->SetSeed(i);}
  double getPropParameter(int ipar) { return pars[ipar];}
  void setCov(covBase *covariance);
  std::string sysType;
#endif
  //saving and reading pars
  void savePars(const char* filename);
  void readPars(const char* filename);
  void printPars();
};

#endif

#ifndef ATMFITPARS_C
#include "atmFitPars.cxx"
#endif


