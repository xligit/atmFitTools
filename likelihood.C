#include "TMath.h"
#include "TH1D.h"

#include <iostream>
#include <math.h>

using namespace std;

/////////////////////////////////////////////////////////////
//Log Likelihood evaluation functions

/////////////////////////////////////////
//likelihood from tn186
double evalLnL(double ndata, double nmc){
  if (nmc==0) return 0.;
  if (ndata==0) return nmc;
  double lnL = (nmc-ndata) + ndata*TMath::Log(ndata/nmc);
  return lnL;
}

//////////////////////////////////////////////////////////
//Use full approximation for comparing to numeric method
double evalLnLRamanujan(double ndata, double nmc){
  if (nmc==0) return 0.;
  if (ndata==0) return nmc;
  double lnL = (nmc-ndata) + ndata*TMath::Log(ndata/nmc)
               + (0.166666)*TMath::Log(ndata*(1+ ((4*ndata)*(1+(2*ndata)))))
               + (0.5)*TMath::Log(TMath::Pi());
  return lnL;
}

//////////////////////////////////////////////////////////////////////
//full numeric calculation
double evalLnLNumeric(double ndata, double mcmean, double mcsig, int ntotpts = 100){
  if (mcsig==0) return 0.;
  double I = 0.;
  double rangemin = fmax((mcmean - (4.*mcsig)),0); 
  double rangemax = (mcmean + (4.*mcsig));
  double dx = (rangemax-rangemin)/(double)ntotpts;
  double x = rangemin;
  for (int ipt = 0;ipt<ntotpts;ipt++){
    double value = TMath::Poisson(ndata,x)*TMath::Gaus(x,mcmean,mcsig,kTRUE);
    I+=(value*dx);
    x+=dx;
  }
  return -1*TMath::Log(I);

}

/////////////////////////////////////////////////////////////////////////
//Choose between numeric and approx to increase speed
double evalLnLFast(double ndata, double mcmean, double mcsig, int ntotpts = 100){
  if ((mcsig/mcmean)<0.05) return evalLnLRamanujan(ndata,mcmean);
  else{
    return evalLnLNumeric(ndata,mcmean,mcsig);
  }
}

void plotL2(double nmc = 4, double errfrac = 0.1, int npts=100){
  double rangemin = 0;
  double rangemax = 10;
  static TH1D* hreal = new TH1D("hreal","hreal",npts,rangemin,rangemax);
  static TH1D* happrox = new TH1D("happrox","happrox",npts,rangemin,rangemax);
  for (int i=0;i<=npts;i++){
    hreal->SetBinContent(i,evalLnLNumeric(hreal->GetBinCenter(i),nmc,errfrac));
    happrox->SetBinContent(i,evalLnLRamanujan(happrox->GetBinCenter(i),nmc));
  }
  hreal->SetLineColor(kBlue);
  hreal->Draw();
  happrox->Draw("same");
  return;
}



void plotL(double ndata = 4, double nmc = 4, int npts=100){
  double rangemin = 0.01;
  double rangemax = 0.8;
  static TH1D* hreal = new TH1D("hreal","hreal",npts,rangemin,rangemax);
  static TH1D* happrox = new TH1D("happrox","happrox",npts,rangemin,rangemax);
  for (int i=0;i<=npts;i++){
    hreal->SetBinContent(i,evalLnLNumeric(ndata,nmc,hreal->GetBinCenter(i)));
    happrox->SetBinContent(i,evalLnLRamanujan(ndata,nmc));
  }
  hreal->SetLineColor(kBlue);
  hreal->Draw();
  happrox->Draw("same");
  return;
}


