#include "TH1F.h"
#include "TRandom2.h"
#include <iostream>
#include <TMath.h>
#include "TString.h"

#define MAXPARS 100

using namespace std;



class anaL{
  void anaL(const char* name);
  TString nameTag;
  void setFCN(void (*fcn)(int &, double* ,double &f,double* ,int));
  void *(Lfcn)(int &, double* ,double &f,double* ,int);
  void profile(int ipar, int npts, float range,const char* hname="");
  //fit imputs
  int Nbin;
  int Ncomp;
  int Natt;
  void setNbin(int nbin) {Nbin=nbin;}
  void setNcomp(int ncomp) {Ncomp=ncomp;}
  void setNatt(int natt) {Natt=natt;}
  //function inputs
  int NPars;
  double* Gout;
  double  Result;
  float   bestPar[10][10][10];
  void    setBestPar(int ibin,int icomp, int iatt, float value);
  double  Par[MAXPARS];
  int     Flg;
};

void setBestPar(int ibin,int icomp, int iatt, float value){
  bestPar[ibin][icomp][iatt] = value;
  return;
}

void anaL:profile(int ipar, int npts, float range, float parval, const char* hname){
  //fix initial parameters
  

  float dx = range/(float)npts;
  float midpoint = range/2.;
  float xx = parval-range;
}

void anaL::setFCN(void (*fcn)(int &, double* ,double &f,double* ,int)){
  Lfcn = fcn;
  return;
}

void anaL::anaL(const char* name){
  nameTag = name;
  cout<<"created likelihood analyzer: "<<nameTag.Data()<<endl;
  return;
}


