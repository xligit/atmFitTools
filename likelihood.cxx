#ifndef LIKELIHOOD_CXX
#define LIKELIHOOD_CXX
#include "TMath.h"

#include "TH1D.h"
#include "TGraph.h"

#include <iostream>
#include <math.h>

using namespace std;

//////////////////////////////////////////////////////////
//Scaled Poisson function since it's not in ROOT
double PoissonS(double nobs, double mean, double scale){
  return TMath::Poisson(nobs/scale,mean);
}

///////////////////////////////////////////////////
//plot scaled Poisson
TGraph* plotPossonS(double mean, double scale, int npts=50){
  double X[npts];
  double Y[npts];
  double rangemin = 0;
  double rangemax = 20;
  double dx = 20/(double)npts;
  double x = 0.;
  for (int ipt=0;ipt<npts;ipt++){
    X[ipt] = x;
    Y[ipt] = PoissonS(x,mean,scale);
    x+=dx;
  }
  TGraph* g = new TGraph(npts,X,Y);
  g->Draw("a*");
  return g;
}


/////////////////////////////////////////////////////////////
//Log Likelihood evaluation functions

/////////////////////////////////////////
//likelihood from tn186
double evalLnL(double ndata, double nmc, double norm = 1.){
  if (nmc==0) return 0.;
  if (ndata==0) return nmc*norm;
  double Nmc = nmc*norm;
  double ngLnL = (Nmc-ndata) + ndata*TMath::Log(ndata/Nmc);
//  if (isnan(ngLnL)){
//    cout<<"WARNING: NAN IN LIKELIHOOOD"<<endl;
//    cout<<"NDATA: "<<ndata<<" NMC: "<<nmc<<" NORM: "<<norm<<endl;
//  }
  return ngLnL;
}


/////////////////////////////////////////
//likelihood from tn186
double evalLnLMean(double ndata, double nmc, double mcerr, double norm = 1.){
  if (nmc==0) return 0.;
  if (ndata==0) return nmc*norm;
  double Nmc = nmc*norm;
  double ss = (mcerr*mcerr/nmc);
  ss = 1.;
  double mean = ((ndata*ss) + (nmc/norm))/((1./norm)+ss);
  double diff = ((1./norm)*ndata) - nmc;
  cout<<"mean: "<<mean<<endl;
  cout<<"ss: "<<ss<<endl;
  double ngLnL = (mean-ndata) + ndata*TMath::Log(ndata/mean);
//  double ngLnL = (diff*diff)/mean;
  return ngLnL;
}



//////////////////////////////////////////////////////////
//Use full approximation for comparing to numeric method
double evalLnLRamanujan(double ndata, double nmc, double norm=1.){
  nmc*=norm;
  if (nmc==0) return 0.;
  if (ndata==0) return nmc*norm;
  double Nmc = nmc*norm;
  double lnL = (Nmc-ndata) + ndata*TMath::Log(ndata/Nmc)
               + (0.166666)*TMath::Log(ndata*(1+ ((4*ndata)*(1+(2*ndata)))))
               + (0.5)*TMath::Log(TMath::Pi());
  return lnL;
}



//////////////////////////////////////////////////////////////////////
//full numeric calculation
double evalLnLNumeric(double ndata, double mcmean, double mcsig, double norm=1., int ntotpts = 100){
//  mcmean*=norm;
//  mcsig*=norm;
  if (mcsig==0) return 0.;
  if (mcmean<0) return 0.;
  double I = 0.;
  double rangemin = fmax((mcmean - (4.*mcsig)),0)*norm; 
  double rangemax = (mcmean + (4.*mcsig))*norm;
  double dx = (rangemax-rangemin)/(double)ntotpts;
  double x = rangemin;
  for (int ipt = 0;ipt<ntotpts;ipt++){
    double value = TMath::Poisson(ndata,x)*TMath::Gaus(x,mcmean*norm,mcsig*norm,kTRUE);
    I+=(value*dx);
    x+=dx;
  }
  return -1*2*TMath::Log(I);

}

/////////////////////////////////////////////////////////////////////////
//
double evalLnLScaled(double ndata, double mcmean, double mcsig,double norm=1.,double ascale=1.){
  if (mcmean==0) return 0.;
 
 // double scale = mcmean/(mcsig*mcsig);
  double scale = (mcsig);

  double ngLnL = evalLnL(ndata,mcmean,norm);
//  cout<<"scale: "<<scale<<endl;
//  ngLnL+= +0.5*TMath::Log(scale)*norm*ascale;
  ngLnL+= 0.5*TMath::Log(scale)*norm*ascale;

  return ngLnL;
}

/////////////////////////////////////////////////////////////////////////
//Assume Gaussian errors
//double evalLnLGauss(double ndata, double mcmean, double mcsig, int ntotpts = 100){
//    return (ndata-mcmean)*(ndata-mcmean)/(2.*mcsig*mcsig);
//}


/////////////////////////////////////////////////////////////////////////
// factor in MC errors
double evalGausChi2WithError(double ndata, double mcmean, double mcsig){

  double sigmasq = ndata + (mcsig*mcsig);
  double diff = (ndata-mcmean);
  return (diff*diff)/(2.*sigmasq);

}

/////////////////////////////////////////////////////////////////////////
//Assume scaled Gaussian errors, and mcmean has been normalized
double evalLnLMyChi2(double ndata, double mcmean, double mcsig, double norm=1.){
    double ss = mcsig;
    double deltasq = (1./(norm*ss))*(mcmean*ss + ndata);
    double diff =(ndata/norm) - mcmean;
    return (diff*diff)/(2.*deltasq);
}

////////////////////////////////////////////////////////////////////////////
double evalLnLDiff(double ndata, double mcmean,  double norm){
  double diff = TMath::Abs(ndata-(mcmean*norm));
  return diff;
}
/////////////////////////////////////////////////////////////////////////
//Assume scaled Gaussian errors, and mcmean has been normalized
double evalLnLMyChi2NoWeight(double ndata, double mcmean, double norm){
    double diff = ndata - (norm*mcmean);
    double deltasq = ( ndata + (norm*mcmean) );
    return -2*TMath::Log(TMath::Gaus(diff,0.,TMath::Sqrt(deltasq),kTRUE));
}



/////////////////////////////////////////////////////////////////////////
//Assume scaled Gaussian errors, and mcmean has been normalized
double evalLnLChi2N(double ndata, double mcmean, double norm){
    double diff = ndata - (mcmean*norm);
    double deltasq = ( ndata );
    return -2*TMath::Log(TMath::Gaus(diff,0.,TMath::Sqrt(deltasq),kTRUE));
}


/////////////////////////////////////////////////////////////////////////
//Assume scaled Gaussian errors
double evalLnLChi2Numeric(double ndata, double mcmean, double mcsig,double norm, int npts = 100){
    double diff = ndata - (norm*mcmean);
    double deltasq = ( (ndata) + (norm*norm*mcsig*mcsig) );
    double sigmasq = ndata;
    double rangemin = fmax(0, mcmean - (5*mcsig))*norm;
    double rangemax = (mcmean + (5*mcsig))*norm;
    double dx = (rangemax-rangemin)/(double)npts;
    double I = 0.;
    double x = rangemin;
    for (int i=0;i<npts;i++){
       double value = TMath::Gaus((x-ndata),0.,deltasq,kTRUE)*TMath::Gaus((x-mcmean),norm*mcmean,mcsig*norm,kTRUE);
       I+=(dx*value);
       x+=dx;
    }
    return -2*TMath::Log(I);
}


void plotLnL(int itype, double nmc = 16, double mcsig = 4., double norm=1, int npts=100){
  double rangemin = 2;
  double rangemax = 40;
  double dx = (rangemax-rangemin)/(double)npts;
  double xx = rangemin;
  static TH1D* hreal = new TH1D("hreal","hreal",npts,rangemin,rangemax);
  double lnL=0.;
  for (int i=0;i<=npts;i++){      
    if (itype==0) lnL=evalLnLRamanujan(xx,nmc,norm); //< tn186
    if (itype==1) lnL=evalLnLNumeric(xx,nmc,mcsig,norm); //< my numerical method
    if (itype==2) lnL=evalLnLMyChi2(xx,nmc,mcsig,norm); //< chi-2 style errors
    if (itype==3) lnL=evalLnLChi2Numeric(xx,nmc,mcsig,norm);
    hreal->SetBinContent(i,lnL);
    xx+=dx;
  }
  hreal->SetLineColor(kBlue);
  hreal->Draw();
  return;
}



void plotL2(double nmc = 4, double errfrac = 0.1, int npts=100){
  double rangemin = 0;
  double rangemax = 20;
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

#endif
