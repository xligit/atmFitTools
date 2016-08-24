#ifndef HISTOTRANSFORMS_C
#define HISTOTRANSFORMS_C

#include "TH1D.h"
#include "TRandom2.h"
#include "TMath.h"
#include "time.h"
#include <iostream>
#include "TCanvas.h"

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
TRandom2* randy = new TRandom2();
#endif


using namespace std;


/////////////////////////////////////////////////////
//returns a parameter for S/N ratio for a histogram
double getNoiseFactor(TH1D* hh){
//  if (hh->GetEntries()==0) return 0.;
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
  if (B==0) return 0;
  S/=B;
  return 1./TMath::Sqrt(S);
}




///////////////////////////////////////////////////////////////////////////////////
//Custom smoothing method
void mySmooth(TH1D* hh,double factor=3.0){

  //////////////////////////////
  //get adjecent bin weights
  
  // use Gaussian weights
  double sigma = factor*getNoiseFactor(hh); //< set sigma using noise
  if (sigma==0) return; //no events in histogram
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

///////////////////////////////////
//Integral of box function
double B(double x,double a, double b){
  if (x<=a) return 0.;
  if (x>b) return 1.;
 // if (a==b) return 0;
  return (x-a)/(b-a);
}


//////////////////////////////////////////////////////////////////////////////////////
//useful test functions
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
///////////////////////////////////////////////////////////////////////////////////////


/*TH1D* convolveHisto(TH1D* hh, double sig, double bias, const char* name=""){
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

*/


/*
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
//hsmear;
  return hsmooth;
}
*/ 

/*
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
*/

/*

//////////////////////////////////////////////////////////////
//performs histogram modifications
void smearThisHisto(TH1D &hh, double spread, double bias=0.){

  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //make temporary clone of input histogram 
  TH1D* htmp = (TH1D*)hh.cxxlone("htmp");

  //apply custom smooth function if statistics are low
  if (hh.GetEntries()<10000.) mySmooth(htmp);

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
//  double mean = hh.GetMean() + (binw/2.);
  double mean = hh.GetMean();
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  double sumw2;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){

 //   if ((htmp->GetBinContent(newbin)<=0)||(htmp->GetBinContent(newbin-1)<=0)||(htmp->GetBinContent(newbin+1)<=0)){
 //     continue;
 //   }

    ////////////////////
    //set weight parameters
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = htmp->GetBinLowEdge(newbin);
 //  ymin = ((binedge-bias)*smear) - shift;
 //   ymin = binedge*smear - shift - bias;
 //   ymax = ymin + (binw*smear);
    ymin  = binedge - bias;
    ymax  = ymin + binw;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = htmp->GetBinContent(oldbin);
      sum+=(weight*htmp->GetBinContent(oldbin));
      binerr += weight*weight*(htmp->GetBinContent(oldbin)*htmp->GetBinContent(oldbin));
      sumw += weight;
    }
    cout<<"sum: "<<sum<<endl;
    cout<<"sumw: "<<sumw<<endl;
    if (sumw<=0.0001) continue;
//    hh.SetBinContent(newbin,sum/sumw);
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(sum));
    double scale = htmp->Integral()/hh.Integral();
    hh.Scale(scale);
  }
  htmp->Delete();
  return;
}

*/


//////////////////////////////////////////////////////////////
//performs histogram modifications
void smearThisHisto(TH1D &hh, double spread, double bias=0.){

  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //make temporary clone of input histogram 
  TH1D* htmp = (TH1D*)hh.Clone("htmp");

  //apply custom smooth function if statistics are low
//  if (hh.GetEntries()<10000.) mySmooth(htmp);

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
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


  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){

 //   if ((htmp->GetBinContent(newbin)<=0)||(htmp->GetBinContent(newbin-1)<=0)||(htmp->GetBinContent(newbin+1)<=0)){
 //     continue;
 //   }

    ////////////////////
    //set weight parameters
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = htmp->GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = htmp->GetBinContent(oldbin);
      sum+=(weight*htmp->GetBinContent(oldbin));
      binerr += weight*weight*(htmp->GetBinContent(oldbin)*htmp->GetBinContent(oldbin));
      sumw += weight;
    }
  //  if (sumw<=0.0001) continue;
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(sum));
//    double scale = htmp->Integral()/hh.Integral();
//    hh.Scale(scale);
  }

  /*
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    ymin = ((htmp->GetBinLowEdge(newbin)-bias)*smear) - shift;
 //   ymin = ((htmp->GetBinLowEdge(newbin))*smear) -bias - shift;
 //   ymax = ((htmp->GetBinLowEdge(newbin)+htmp->GetBinWidth(newbin)-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+htmp->GetBinWidth(oldbin));
      weight = B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      sum+=(weight*htmp->GetBinContent(oldbin));
    }
    hh.SetBinContent(newbin,sum);
  }
  */

  double scale = htmp->Integral()/hh.Integral();
  hh.Scale(scale);

  htmp->Delete();
  return;
}



//////////////////////////////////////////////////////////////
//Test function to compare to event-by-event
void compareEvtByEvt(int nevts, double mean0, double sig0,double scale, double bias){

  //nominal histogram
  TH1D* h0 = testBumpD(nevts, sig0, mean0, "nominal");
  double mean = h0->GetMean();

  //transformed histogram
  TH1D* htransform =  testBumpD(nevts, sig0, mean0, "transformed");
  htransform->SetLineColor(kRed);
  smearThisHisto(*htransform,scale,bias);

  //event by event histogram
  TH1D* hebe=  testBumpD(nevts, sig0, mean0, "evntbyevt");
  hebe->Reset();
  hebe->SetLineColor(kBlue);
  for (int i=0;i<nevts;i++){
    double xx = randy->Gaus(mean0,sig0); //< get random throw
    xx *= scale;
    xx += bias;
   // xx += bias*scale;
    xx -= scale*mean;
    xx += mean;
    hebe->Fill(xx);
  }

  TCanvas* cc = new TCanvas("cc","cc",800,700);
  h0->Draw();
  htransform->Draw("same");
  hebe->Draw("same");
  cc->Print("compare.png"); 

  return;

}


//smear it faster without mean
void smearThisHistoFast(TH1D &hh, double* hcontent, double spread,  double bias, double normscale=1.){

  int nbins=hh.GetNbinsX();
 // double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double smear = 1./spread;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
//    ymin = binedge - bias;
//    ymax = ymin + binw;
    ymin = ((binedge-bias)*smear);
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*binc;
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(binerr));
  }
  return;
}



//smear it faster without mean
void smearThisHistoFastBias(TH1D &hh, double* hcontent, double bias, double normscale=1.){

  int nbins=hh.GetNbinsX();
 // double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
    ymin = binedge - bias;
    ymax = ymin + binw;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*binc;
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(binerr));
  }
  return;
}



//smear it faster
void smearThisHistoFastMean(TH1D &hh, double* hcontent, double spread, double mean, double bias, double normscale=1.){

  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //make temporary clone of input histogram 
//  TH1D* htmp = (TH1D*)hh.Clone("htmp");

  //apply custom smooth function if statistics are low
  // if (hh.GetEntries()<10000.) mySmooth(h);

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double binwidth;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
//  double mean = hh.GetMean() + (binw/2.);
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  double sumw2;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
//    ymin = (binedge-bias)*
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*(binc*binc);
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum);
  }
  double newintegral = hh.Integral();
  double scale = oldintegral/newintegral;
  hh.Scale(scale*normscale);
//  mySmooth(&hh);
  return;
}

TH1D* testsmearmean(double bias, double smear, int nev, const char* name = "testsmear.png"){
 
  double mean = 2.0;

  //canvas
  TCanvas* cc = new TCanvas("cc","cc",800,700);
  

  //make a test histo
  TH1D* hh = testBumpD(nev,1.0,2.0); 
  hh->Sumw2();
  TH1D* horig = (TH1D*)hh->Clone("orig");

  horig->SetLineColor(kRed);
  horig->Draw("h");

  //make content
  double hcont[1000]; 
  for (int i=0; i<hh->GetNbinsX();i++){
    hcont[i] = hh->GetBinContent(i);
  }

  //smear it
  smearThisHistoFastMean(*hh, hcont,smear,mean, bias,1);

  hh->Draw("sameh");

  cc->Print(name);

  return hh;
  

}


TH1D* testsmear(double smear, double bias, int nev, const char* name = "testsmear.png"){

  //canvas
  TCanvas* cc = new TCanvas("cc","cc",800,700);
  

  //make a test histo
  TH1D* hh = testBumpD(nev,1.0,1.0); 
  hh->Sumw2();
  TH1D* horig = (TH1D*)hh->Clone("orig");

  horig->SetLineColor(kRed);
  horig->Draw("h");

  //make content
  double hcont[1000]; 
  for (int i=0; i<hh->GetNbinsX();i++){
    hcont[i] = hh->GetBinContent(i);
  }

  //smear it
  smearThisHistoFast(*hh, hcont,smear, bias, 1);

  hh->Draw("sameh");

  cc->Print(name);

  return hh;
  

}


/*
//smear it faster
void smearThisHistoFast(TH1D &hh, double* hcontent, double spread, double bias=0.){
  if (hh.Integral()<1e-4) return;
  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //make temporary clone of input histogram 
//  TH1D* htmp = (TH1D*)hh.Clone("htmp");

  //apply custom smooth function if statistics are low
  // if (hh.GetEntries()<10000.) mySmooth(htmp);

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
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

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*(binc*binc);
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum/sumw);
  }
  double newintegral = hh.Integral();
  double scale = oldintegral/newintegral;
  hh.Scale(scale);
  return;
}
*/

/*
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

*/


double testtime(){
  int ntry = 25000;
  clock_t t1,t2;
  TH1D* hb = testBump(1000);
  TH1D* hmod = (TH1D*)hb->Clone("htmp");
  t1=clock();
  for (int i=0;i<ntry;i++){
    smearThisHisto(*hb,1.1,1.2);
  }
  t2=clock();
  double diff = ((double)t2-(double)t1)/((double)ntry);
  cout<<"time: "<<diff<<endl;
  return diff;
  
}




#endif

