#ifndef ATMFITPARS_C
#define ATMFITPARS_C

#include "shared.h"

using namespace std;

//class containing all atmospheric fit pars
class atmFitPars{
  public:
  
  atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst=0);
 
  //numbers of various parametrs
  int nSamples;
  int nBins;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nTotPars;

  //parameter values
  float histoPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  float histoParUncLo[NBINMAX][NCOMPMAX][NATTMAX][2];
  float histoParUncHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  float sysPar[NSYSPARMAX];
  float sysParUnc[NSYSPARMAX];
  float pars[4000];
  int   fixPar[4000];
  float bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2];
  //normalization
  float norm;  
  void setNorm(float x){norm=x;}
  void initPars();
  int getParIndex(int ibin, int icomp, int iatt, int imod){return parIndex[ibin][icomp][iatt][imod];}
  void setParameter(int ipar, float value);
  void setSysParameter(int ipar, float value);
  void setParameter(int ibin, int icomp, int iatt, int imod, float value); 
  void setSysParUnc(int isys,float value){sysParUnc[isys]=value;}
  void fixParameter(int ipar);
  void fixParameter(int ibin,int icomp,int iatt, int imod);
  int  checkFixFlg(int ibin,int icomp,int iatt, int imod);
  void printParValues();
  int binOfPar[4000];
  int compOfPar[4000];
  int attOfPar[4000];
  int typeOfPar[4000];
};

void atmFitPars::printParValues(){
  for (int ipar=0;ipar<nTotPars;ipar++){
    cout<<"par "<<ipar<<": "<<pars[ipar]<<endl;
  }
}

int atmFitPars::checkFixFlg(int ibin, int icomp, int iatt, int imod){
  return fixPar[parIndex[ibin][icomp][iatt][imod]];
}

void atmFitPars::fixParameter(int ibin,int icomp,int iatt, int imod){
  fixPar[getParIndex(ibin,icomp,iatt,imod)] = 1;
  return;
}

void atmFitPars::fixParameter(int ipar){
  fixPar[ipar] = 1;
  return;
}

void atmFitPars::setSysParameter(int ipar, float value){
  sysPar[ipar]=value;
  pars[nTotPars-nSysPars+ipar]=value;
  return;
}

void atmFitPars::setParameter(int ibin, int icomp, int iatt, int itype, float value){
  histoPar[ibin][icomp][iatt][itype] = value;
  pars[parIndex[ibin][icomp][iatt][itype]] = value;
  return;
}

void atmFitPars::setParameter(int ipar, float value){
  pars[ipar]=value;
  if (ipar>=(nTotPars-nSysPars)) sysPar[nTotPars-ipar-1] = value;
  else{
    histoPar[binOfPar[ipar]][compOfPar[ipar]][attOfPar[ipar]][typeOfPar[ipar]]=value;
  }
  return;
}

atmFitPars::atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst){
  nSamples = isamp;
  nBins    = ibin;
  nComponents    = icomp;
  nAttributes   = iatt;
  nSysPars = nsyst;
  initPars();
}

void atmFitPars::initPars(){
  int index = 0;
  for (int ibin=0;ibin<nBins;ibin++){
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        histoPar[ibin][icomp][iatt][0]=1.0;
        pars[index]=1.0;
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=0;
        fixPar[index]=0;
        parIndex[ibin][icomp][iatt][0]=index;  
        index++;
        histoPar[ibin][icomp][iatt][1]=0.0;
        parIndex[ibin][icomp][iatt][1]=index;
        pars[index]=0.0;
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=1;
        fixPar[index]=0;
        index++;
      }
    }
  }
  for (int isys=0;isys<nSysPars;isys++){
    sysPar[isys]=1.;
    pars[index]=1.;
    index++;
  }
  nTotPars = index;
  cout<<"Total number of fit parameters: "<<nTotPars<<endl;
  for (int kpar=0;kpar<nTotPars;kpar++){
    cout<<"par "<<kpar<<" value: "<<pars[kpar]<<endl;
  }
  return;
}



          
#endif
