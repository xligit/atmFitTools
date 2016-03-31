#include "TH1D.h"

using namespace std;

//Integral of box function
double B(double x,double a, double b){
  if (x<=a) return 0.;
  if (x>=b) return 1.;
  double slope = 1./(b-a);
  return (x-a)/(b-a);
}


//Returns a smeared histogram
//h -> Initial histo
//spread -> Factor to determine width smearing. <1 will shrink histo, >1 will stretch
//bias -> Adds bias to histogram 
TH1D* smearIt(TH1D* h,double spread, double bias=0.){
  TH1D* hsmear = (TH1D*)h->Clone("hsmear");
  hsmear->Reset();
  int nbins=hsmear->GetNbinsX();
  double sum;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
  double mean = h->GetMean();
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
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


