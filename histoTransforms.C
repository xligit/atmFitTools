#ifndef HISTOTRANSFORMS_C
#define HISTOTRANSFORMS_C

#include "TH1F.h"
#include "TRandom2.h"
#include "TMath.h"
#include "time.h"
#include <iostream>

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
TRandom2* randy = new TRandom2();
#endif


using namespace std;

//Integral of box function
float B(float x,float a, float b){
  if (x<=a) return 0.;
  if (x>=b) return 1.;
  return (x-a)/(b-a);
}

TH1F* testBump(int nev,double sig=1.0,double mean=0.0){
  TH1F* h = new TH1F("testbump","testbump",100,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1F* convolveHisto(TH1F* hh, float sig, float bias, const char* name=""){
  //name setup
  TString hname = "convhist";
  hname.Append(name);
 
  //make convolved histogram;
  TH1F* hconv;
  hconv = (TH1F*)hh->Clone(hname.Data());
  hconv->Reset();

  //vars for calculation
  int nbinsx = hh->GetNbinsX();
  float xx;
  float area;
  
  float bmin;
  float bmax;
  float gweight; //integral of gaussian in bin
  for (int ibin=1;ibin<=nbinsx;ibin++){
    //calculate area overlap for each bin:
    xx = hh->GetBinCenter(ibin);
    area = 0.;
    for (int jbin=1;jbin<=nbinsx;jbin++){
      bmin = hh->GetBinLowEdge(jbin);
      bmax =  bmin+hh->GetBinWidth(jbin);
      gweight = 1.;
      gweight = TMath::Erf((bmax -(xx-bias))/(sig*sqrt(2.))) - TMath::Erf((bmin
-(xx-bias))/(sig*sqrt(2.)));
      area+=gweight*hh->GetBinContent(jbin);
    }
    hconv->SetBinContent(ibin,area); 
  }

  //normalize
  float norm = (float)hh->Integral()/(float)hconv->Integral();
  hconv->Scale(norm);
//  hconv->SetLineColor(kRed);
//  hconv->Draw();
//  hh->Draw("same");
  return hconv;
}

//Returns a smeared histogram
//h -> Initial histo
//spread -> Factor to determine width smearing. <1 will shrink histo, >1 will stretch
//bias -> Adds bias to histogram 


TH1F* smearIt(TH1F* h,float spread, float bias=0.){
  TH1F* hsmear = (TH1F*)h->Clone("hsmear");
  hsmear->Reset();
  int nbins=hsmear->GetNbinsX();
  float sum;
  float weight;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  if (spread==0) return hsmear;
  float smear = 1./spread;
  float mean = h->GetMean() + (h->GetBinWidth(1)/2.);
  float shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    ymin = ((hsmear->GetBinLowEdge(newbin)-bias)*smear) - shift;
    ymax = ((hsmear->GetBinLowEdge(newbin)+hsmear->GetBinWidth(newbin)-bias)*smear) - shift;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = h->GetBinLowEdge(oldbin);
      xmax = (xmin+h->GetBinWidth(oldbin));
      weight = B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      sum+=(weight*h->GetBinContent(oldbin));
    }
    hsmear->SetBinContent(newbin,sum);
  }
  if (hsmear->Integral()>0.) hsmear->Scale(h->Integral()/hsmear->Integral());
  hsmear->SetLineColor(kMagenta);
 // TH1F* hsmooth = convolveHisto(hsmear,(hsmear->GetBinWidth(1)/2.),0.);
  /* // for debugging
  if (spread<=1.0){
    hsmear->Draw();
    h->Draw("same");
  }
  else{
    h->Draw();
    hsmear->Draw("same");
  }
  */
  return hsmear;
//  return hsmooth;
}



void smearThisHisto(TH1F &hh, float spread, float bias=0.){
  if (spread==0) return;
  TH1F* htmp = (TH1F*)hh.Clone("htmp");
  int nbins=hh.GetNbinsX();
  float binw = hh.GetBinWidth(1);
  float binedge;
  float sum;
  float weight;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float smear = 1./spread;
  float mean = hh.GetMean() + (binw/2.);
  float shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    binedge = htmp->GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    //ymax = ((binedge+binw-bias)*smear) - shift;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight = B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      sum+=(weight*htmp->GetBinContent(oldbin));
    }
    hh.SetBinContent(newbin,sum);
  }
  if (hh.Integral()>0.) hh.Scale(htmp->Integral()/hh.Integral());
  htmp->Delete();
  return;
}


void convolveThisHisto(TH1F &hh, float sig, float bias){
 
  //make convolved histogram;
  TH1F* htmp;
  htmp = (TH1F*)hh.Clone("tmphisto");

  //vars for calculation
  int nbinsx = hh.GetNbinsX();
  float xx;
  float area;
  float bmin;
  float bmax;
  float binw = htmp->GetBinWidth(1);
  float sqrt2 = sqrt(2.);
  float gweight; //integral of gaussian in bin
  float rms0 = (float)hh.GetRMS();
  for (int ibin=1;ibin<=nbinsx;ibin++){
    //calculate area overlap for each bin:
    xx = htmp->GetBinCenter(ibin);
    area = 0.;
    for (int jbin=1;jbin<=nbinsx;jbin++){
      bmin = htmp->GetBinLowEdge(jbin);
      bmax =  bmin+binw;
      gweight = 1.;
      gweight = TMath::Erf((bmax -(xx-bias))/(sig*sqrt2)) - TMath::Erf((bmin-(xx-bias))/(sig*sqrt2));
      area+=gweight*htmp->GetBinContent(jbin);
    }
    hh.SetBinContent(ibin,area); 
  }

  //normalize
  float norm = (float)htmp->Integral()/(float)hh.Integral();
  hh.Scale(norm);
  float rms1 = (float)hh.GetRMS();
  float correction = rms0/(rms1);
 // cout<<"smear factor: "<<correction<<endl;
  smearThisHisto(hh,correction);
//  hconv->SetLineColor(kRed);
//  hconv->Draw();
//  hh->Draw("same");
  htmp->Delete();
  return;
}







void smearHisto(TH1F &hi,TH1F &hf,float spread, float bias=0.){
  if (spread==0) return;
  int nbins=hi.GetNbinsX();
  float binw = hi.GetBinWidth(1);
  float binedge;
  float sum;
  float weight;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float smear = 1./spread;
  float mean = hi.GetMean() + (binw/2.);
  float shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    binedge = hi.GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    //ymax = ((binedge+binw-bias)*smear) - shift;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hi.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight = B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
    //  sum+=(weight*hi.GetBinContent(oldbin));
      sum+=(weight*(hi.GetBinContent(oldbin)+((0.5)*(hi.GetBinContent(oldbin-1))+(0.5*hi.GetBinContent(oldbin+1)))));
    }
    hf.SetBinContent(newbin,sum);
  }
  if (hf.Integral()>0.) hf.Scale(hi.Integral()/hf.Integral());
  return;
}




float testtime(){
  int ntry = 25000;
  clock_t t1,t2;
  TH1F* hb = testBump(1000);
  TH1F* hmod = (TH1F*)hb->Clone("htmp");
  t1=clock();
  for (int i=0;i<ntry;i++){
    smearHisto(*hb,*hmod,1.1,1.2);
  }
  t2=clock();
  float diff = ((float)t2-(float)t1)/((float)ntry);
  cout<<"time: "<<diff<<endl;
  return diff;
  
}



#endif

