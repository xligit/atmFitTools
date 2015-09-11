#include "TH1F.h"

using namespace std;

//Integral of box function
float B(float x,float a, float b){
  if (x<=a) return 0.;
  if (x>=b) return 1.;
  float slope = 1./(b-a);
  return (x-a)/(b-a);
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
  float smear = 1./spread;
  float mean = h->GetMean();
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
  float binx;
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




