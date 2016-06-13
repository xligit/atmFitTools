#ifndef ATMFITPARS_H
#define ATMFITPARS_H

#include "shared.h"
#include "TRandom2.h"
#include "sharedPars.cxx"

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

  ///////////////////////////////////////////////////////////////
  //numbers of various parametrs
  int nSamples;
  int nBins;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nTotPars;
  int flgUseNormPars;
  sharedPars* runpars;

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

  //saving and reading pars
  void savePars(const char* filename);
  void readPars(const char* filename);
  void printPars();
};

#endif

#ifndef ATMFITPARS_C
#include "atmFitPars.cxx"
#endif


