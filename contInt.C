#include "TH1F.h"
#include <cmath>
#include <iostream>

using namespace std;


//smoothly integrate histogram assuming lines between bins


float fullBinIntegral(TH1F &hh, int ibin){
  float binhw = hh.GetBinWidth(ibin)*0.5; //assume constant bin width
  float slopei = (hh.GetBinContent(ibin)-hh.GetBinContent(ibin-1))/(2.*binhw);
  float slopef = (hh.GetBinContent(ibin+1)-hh.GetBinContent(ibin))/(2.*binhw);
  float integral = 0.;
  float content = hh.GetBinContent(ibin);
  cout<<"slopii: "<<slopei<<endl;
  cout<<"slopif: "<<slopef<<endl; 
  integral+=(binhw*binhw*(slopef-slopei))/2. + (content*(binhw*2.));
  return integral;
}

float partialBinIntegralHigh(TH1F &hh,  float x){
  int ibin = hh.FindBin(x);
  float slopei = (hh.GetBinContent(ibin)-hh.GetBinContent(ibin-1));
  float slopef = (hh.GetBinContent(ibin+1)-hh.GetBinContent(ibin));
  float ww;
  float bincenter = hh.GetBinCenter(ibin);
  float binhw = hh.GetBinWidth(ibin)*0.5; //assume constant bin width
  cout<<"binhw: "<<binhw<<endl;
  float integral = 0.;
  float content = hh.GetBinContent(ibin);
  cout<<"content: "<<content<<endl;  
  //calculate slopes
  cout<<"center: "<<bincenter<<endl;
  ww = fabs(x-bincenter);
  cout<<"ww: "<<ww<<endl;
  slopei/=(binhw*2.);
  slopef/=(binhw*2.);
  cout<<"slopii: "<<slopei<<endl;
  cout<<"slopif: "<<slopef<<endl; 
  if (x<bincenter){
    integral+= (binhw-ww)*(content - (slopei*0.5*(ww+binhw)) );
  }
  else{
    integral+=(content*(ww) + (ww*ww*slopef)/2. );
    integral+=((content*(binhw))-(binhw*binhw*slopei/2.));
  }
  return integral;
}

float partialBinIntegralLow(TH1F &hh, float x){
  int ibin = hh.FindBin(x);
  float slopei = (hh.GetBinContent(ibin)-hh.GetBinContent(ibin-1));
  float slopef = (hh.GetBinContent(ibin+1)-hh.GetBinContent(ibin));
  float ww;
  float bincenter = hh.GetBinCenter(ibin);
  float binhw = hh.GetBinWidth(ibin)*0.5; //assume constant bin width
  cout<<"binhw: "<<binhw<<endl;
  float integral = 0.;
  float content = hh.GetBinContent(ibin);
  cout<<"content: "<<content<<endl;  
  cout<<"center: "<<bincenter<<endl;
  //calculate slopes
  ww = fabs(bincenter-x);
  cout<<"ww: "<<ww<<endl;
  slopei/=(binhw*2.);
  slopef/=(binhw*2.);
  cout<<"slopii: "<<slopei<<endl;
  cout<<"slopif: "<<slopef<<endl; 
  if (x>bincenter){
    integral+= (binhw-ww)*(content+ (slopef*0.5*(ww+binhw)) );
  }
  else{
    integral+=(content*(ww) - (ww*ww*slopei)/2. );
    integral+=((binhw)*(binhw)*slopef)/2. + (content*(binhw));
  }
  cout<<"low component:"<<integral<<endl;
  return integral;
}



float contInt(TH1F &hh, float x1, float x2, float scale=1., float bias=0.){
  int bin1 = hh.FindBin(x1);
  int bin2 = hh.FindBin(x2);
  float integral = 0.;
  if (bin1!=bin2){
    integral+=partialBinIntegralLow(hh,x1);
    integral+=partialBinIntegralHigh(hh,x2);
  }
  else{
    integral+= fullBinIntegral(hh,bin1)-partialBinIntegralLow(hh,x2)
               -partialBinIntegralHigh(hh,x1);
  }
  for (int ibin=(bin1+1);ibin<(bin2-1);ibin++){
    integral+=fullBinIntegral(hh,ibin);
  }
  return integral;
}



TH1F* drawTrans(TH1F &hh, int npts = 10, float scale=1.,float bias=0.){
  TH1F* hnew = new TH1F("hnew","hnew",npts,hh.GetBinLowEdge(1),(hh.GetBinWidth(1)+hh.GetBinLowEdge(hh.GetNbinsX())));
  float xlow;
  float xhigh;
  for (int ibin=1;ibin<=npts;ibin++){
    xlow = hnew->GetBinLowEdge(ibin);
    xhigh = xlow+hnew->GetBinWidth(ibin);
    xlow = (xlow - bias)/scale;
    xhigh = (xhigh - bias)/scale;
    hnew->SetBinContent(ibin,contInt(hh,xlow,xhigh));
  }
  hnew->Scale(hh.Integral()/hnew->Integral());
  hnew->Draw();
  return hnew;
}
