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
  hsmear->Scale(h->Integral()/hsmear->Integral());
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


