#include "splineManager.h"

void getEventWeight(){
  float ww=0.;
  if (fqReader->mode==1)ww*=0.9;
  return ww;
}

void splineManager::splineManager(const char* name, int nsample, int nbin, int ncomp, int natt, int nsyst){
  int nBin=nbin;
  int nComp=ncomp;
  int nAtt=natt;
  int nSyst=nsyst;
}

void splineManager::setSysUnc(int ipar, float val){
  sysUnc[ipar]=val;
  return;
}

void splineManager::setSysPar(int ipar, float val){
  sysPar[ipar]=val
  return;
}

void splineManager::setMCTree(TTree* tr){
  mcTree = tr;
}






