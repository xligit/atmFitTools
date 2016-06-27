#ifndef MASKTOOLS_C
#define MASKTOOLS_C

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include <iostream>
#include "FVCalculators.cxx"
#ifndef T2K
#include "fqEvent.h"
#else
#include "t2kfqEvent.h"
#endif

using namespace std;

int passMask(TH1D* hmask, double value){
  int ibin = hmask->FindBin(value);
  return (int)hmask->GetBinContent(ibin);
}

// a class to mask out spikes in hybrid pi0 data samples
class masktools{
  public:

  // constructor
  masktools();

  // chain with all files
  TChain* chdata;

  // number of bins in histogram
  int nbins;

  // threshold for masking
  double thresh;

  // histgrams
  TH1D* hwall;
  TH1D* hmask;

  // make a mask and save it
  void makethismask(const char* filename);

};

void masktools::makethismask(const char* filename){
  
  // make histograms
  hwall = new TH1D("hwall","hwall",nbins,0,1000);
  hmask = (TH1D*)hwall->Clone("hmask");

  // fill histogram
#ifndef T2K
  fqEvent* fqevent = new fqEvent(chdata);
#else
  t2kfqEvent *fqevent = new t2kfqEvent(chdata);
#endif
  int nev = chdata->GetEntries();
  TVector3* vpos = new TVector3();
  for (int iev=0; iev<nev; iev++){
    chdata->GetEntry(iev);
    vpos->SetXYZ(fqevent->fq1rpos[0][2][0],
                 fqevent->fq1rpos[0][2][1],
                 fqevent->fq1rpos[0][2][2]);
    double wall = calcWall(vpos);
    hwall->Fill(wall);
  }

  // make mask 
  for (int ibin=1; ibin<=hwall->GetNbinsX(); ibin++){
     double content = hwall->GetBinContent(ibin);
     if (content<thresh) hmask->SetBinContent(ibin,1.);
     else{
       hmask->SetBinContent(ibin,0.);
     }
  }

  // save mask
  hmask->SaveAs(filename);

  return;
}

masktools::masktools(){
  
}











#endif
