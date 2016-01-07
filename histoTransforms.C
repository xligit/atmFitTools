#ifndef HISTOTRANSFORMS_C
#define HISTOTRANSFORMS_C

#include "TH1D.h"
#include "TRandom2.h"
#include "TMath.h"
#include "time.h"
#include <iostream>

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
TRandom2* randy = new TRandom2();
#endif


using namespace std;


double getNoiseFactor(TH1D* hh){
  int nbins = hh->GetNbinsX();
  double S = 0.;
  double B = 0;
  for (int ibin = 1;ibin<=nbins;ibin++){
    double content = hh->GetBinContent(ibin);
    if (content>0.){
      S+=content;
      B++;
    }
  }
  S/=B;
  return 1./TMath::Sqrt(S);
}




///////////////////////////////////////////////////////////////////////////////////
//Custom smearing method
void mySmooth(TH1D* hh,double factor=3.){

  //////////////////////////////
  //get adjecent bin weights
  
  // use Gaussian weights
  double sigma = factor*getNoiseFactor(hh); //< set sigma using noise
  double w2 = TMath::Gaus(2,0,sigma,1);
  double w1 = TMath::Gaus(1,0,sigma,1);
  double w0 = TMath::Gaus(0,0,sigma,1);

  ////////////////////////////////
  //clone in original histogram
  TH1D* htmp = (TH1D*)hh->Clone("tmphistosmooth");

  ///////////////////////////////////////////////
  //set new histogram contents from adjecent bins
  double newcontent;
  for (int ibin=0;ibin<hh->GetNbinsX();ibin++){
    newcontent = 0.;
    newcontent+=(htmp->GetBinContent(ibin-2)*w2);
    newcontent+=(htmp->GetBinContent(ibin-1)*w1);
    newcontent+=(htmp->GetBinContent(ibin)*w0);
    newcontent+=(htmp->GetBinContent(ibin+1)*w1);
    newcontent+=(htmp->GetBinContent(ibin+2)*w2);
    hh->SetBinContent(ibin,newcontent);
    hh->SetBinError(ibin,TMath::Sqrt(newcontent));
  }

  //normalize histogram
  double normscale = htmp->Integral()/hh->Integral();
  hh->Scale(normscale); 
  htmp->Delete(); 
  return; 
}

//Integral of box function
double B(double x,double a, double b){
  if (x<=a) return 0.;
  if (x>=b) return 1.;
  return (x-a)/(b-a);
}

TH1D* testBumpD(int nev,double sig=1.0,double mean=0.0,const char* name="testbumb"){
  TH1D* hset;
  hset->SetDefaultSumw2(kTRUE);
  TH1D* h = new TH1D(name,name,40,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1D* testBump(int nev,double sig=1.0,double mean=0.0,const char* name="testbumb"){
  TH1D* hset;
  hset->SetDefaultSumw2(kTRUE);
  TH1D* h = new TH1D(name,name,40,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1D* testTable(int nev,double width=4.0,double mean=0.0){
  TH1D* h = new TH1D("testbump","testbump",50,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill((randy->Rndm()*width) - (width/2.)+mean);
  }
  return h;
}



TH1D* convolveHisto(TH1D* hh, double sig, double bias, const char* name=""){
  //name setup
  TString hname = "convhist";
  hname.Append(name);
 
  //make convolved histogram;
  TH1D* hconv;
  hconv = (TH1D*)hh->Clone(hname.Data());
  hconv->Reset();

  //vars for calculation
  int nbinsx = hh->GetNbinsX();
  double xx;
  double area;
  
  double bmin;
  double bmax;
  double gweight; //integral of gaussian in bin
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
  double norm = (double)hh->Integral()/(double)hconv->Integral();
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
  if (spread==0) return hsmear;
  double smear = 1./spread;
  double mean = h->GetMean() + (h->GetBinWidth(1)/2.);
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
  if (hsmear->Integral()>0.) hsmear->Scale(h->Integral()/hsmear->Integral());
  hsmear->SetLineColor(kMagenta);
 // TH1D* hsmooth = convolveHisto(hsmear,(hsmear->GetBinWidth(1)/2.),0.);
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



void convolveThisHisto(TH1D &hh, double sig, double bias){
 
  //make convolved histogram;
  TH1D* htmp;
  htmp = (TH1D*)hh.Clone("tmphisto");

  //vars for calculation
  int nbinsx = hh.GetNbinsX();
  double xx;
  double area;
  double bmin;
  double bmax;
  double binw = htmp->GetBinWidth(1);
  double sqrt2 = sqrt(2.);
  double gweight; //integral of gaussian in bin
  double rms0 = (double)hh.GetRMS();
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
  double norm = (double)htmp->Integral()/(double)hh.Integral();
  hh.Scale(norm);
  double rms1 = (double)hh.GetRMS();
  double correction = rms0/(rms1);
 // cout<<"smear factor: "<<correction<<endl;
//  smearThisHisto(hh,correction);
//  hconv->SetLineColor(kRed);
//  hconv->Draw();
//  hh->Draw("same");
  htmp->Delete();
  return;
}


void smearThisHisto(TH1D &hh, double spread, double bias=0.){
  if (spread==0){
    cout<<"smearThisHisto: cannont smear with 0 spread parameter"<<endl;
    return;
  }
  //cout<<"cloneing input"<<endl;
  hh.SetDefaultSumw2(kTRUE);
  TH1D* htmp = (TH1D*)hh.Clone("htmp");
  if (hh.GetEntries()<10000.) mySmooth(htmp);
 // float entriesperbin = htmp->GetEntries()/(float)htmp->GetNbinsX();
  //cout<<entriesperbin<<endl;
 // if (entriesperbin<1) htmp->Smooth(1);
 // convolveThisHisto((*htmp),(htmp->GetBinWidth(1)/1.),0.);
  int nbins=hh.GetNbinsX();
  double binw = hh.GetBinWidth(1);
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
  double mean = hh.GetMean() + (binw/2.);
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  double sumw2;
  double sumw;
  for (int newbin=1;newbin<=nbins;newbin++){
    if ((htmp->GetBinContent(newbin)<=0)||(htmp->GetBinContent(newbin-1)<=0)||(htmp->GetBinContent(newbin+1)<=0)){
      continue;
    }
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = htmp->GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    //ymax = ((binedge+binw-bias)*smear) - shift;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax) ;
      double binc = htmp->GetBinContent(oldbin);
      sum+=(weight*htmp->GetBinContent(oldbin));

//      if (binc>100)  sum+=(weight*htmp->GetBinContent(oldbin));
//      else{
//        sum+=(weight*(0.3333333333)*(htmp->GetBinContent(oldbin-1)+htmp->GetBinContent(oldbin)+htmp->GetBinContent(oldbin+1)));
//      }
      binerr += weight*weight*(htmp->GetBinContent(oldbin)*htmp->GetBinContent(oldbin));
      sumw += weight;
    }
    if (sumw<=0) continue;
    hh.SetBinContent(newbin,sum/sumw);
    double ss = binerr/(sum*sum);
    //set bin uncertainty..
    hh.SetBinError(newbin,(sum*sum*sum)/(binerr));
  //  hh.SetBinError(newbin,ss);
  }
  htmp->Delete();
  return;
}


void smearThisHistoD(TH1D &hh, double spread, double bias=0.){
  if (spread==0){
    cout<<"smearThisHisto: cannont smear with 0 spread parameter"<<endl;
    return;
  }
  //cout<<"cloneing input"<<endl;
  hh.SetDefaultSumw2(kTRUE);
  TH1D* htmp = (TH1D*)hh.Clone("htmp");
 //  htmp->Smooth(0);
 // convolveThisHisto((*htmp),(htmp->GetBinWidth(1)/1.),0.);
  int nbins=hh.GetNbinsX();
  double binw = hh.GetBinWidth(1);
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
  double mean = hh.GetMean() + (binw/2.);
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing

  //cout<<"calculating new bin content"<<endl;
  for (int newbin=1;newbin<=nbins;newbin++){
    //////////////////////////////////
    //it is necessary to "mask out" bins where the mc predicts zero events
    //what is the uncertainty on the mean in this case?
//    if (htmp->GetBinContent(newbin)<=0.){
    if ((htmp->GetBinContent(newbin)<=0.001)||(htmp->GetBinContent(newbin-1)<=0.001)||(htmp->GetBinContent(newbin+1)<=0.001)){
//    if (1){
  //    cout<<"no calculation for bin: "<<newbin<<endl;
      hh.SetBinContent(newbin,htmp->GetBinContent(newbin));
      hh.SetBinError(newbin,htmp->GetBinError(newbin));
      continue;
    }
  //  */
    sum = 0.;
    binerr=0.;
    binedge = htmp->GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    //ymax = ((binedge+binw-bias)*smear) - shift;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight = B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      sum+=(weight*htmp->GetBinContent(oldbin));
      //calculate new bin error
   //   binerr += weight*weight*(htmp->GetBinContent(oldbin));
      binerr += weight*weight*(htmp->GetBinError(oldbin))*(htmp->GetBinError(oldbin));
  //    cout<<"content: "<<htmp->GetBinContent(oldbin)<<" error: "<<htmp->GetBinError(oldbin)*htmp->GetBinError(oldbin)<<endl;
     // sum+=(weight*(htmp->GetBinContent(oldbin)+(0.5*htmp->GetBinContent(oldbin-1))+(0.5*htmp->GetBinContent(oldbin+1))));
      
    }
    hh.SetBinContent(newbin,sum);
    //set bin uncertainty..
  //  hh.SetBinError(newbin,TMath::Sqrt(binerr));
    hh.SetBinError(newbin,TMath::Sqrt(sum));

  }
  if (hh.Integral()>0.) hh.Scale(htmp->Integral()/hh.Integral());
//  hh.Smooth(5);
  htmp->Delete();
  return;
}




void smearHisto(TH1D &hi,TH1D &hf,double spread, double bias=0.){
  if (spread==0) return;
  int nbins=hi.GetNbinsX();
  double binw = hi.GetBinWidth(1);
  double binedge;
  double sum;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
  double mean = hi.GetMean() + (binw/2.);
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
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




double testtime(){
  int ntry = 25000;
  clock_t t1,t2;
  TH1D* hb = testBump(1000);
  TH1D* hmod = (TH1D*)hb->Clone("htmp");
  t1=clock();
  for (int i=0;i<ntry;i++){
    smearHisto(*hb,*hmod,1.1,1.2);
  }
  t2=clock();
  double diff = ((double)t2-(double)t1)/((double)ntry);
  cout<<"time: "<<diff<<endl;
  return diff;
  
}




#endif

