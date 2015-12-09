#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#include "TFile.h"

#define NMCMCPARS 100

using namespace std;

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
TRandom2* randy = new TRandom2();
#endif

//class to manage a Markov Chain Monte Carlo
class markovTools{
   public:
   //constructor
   markovTools(int npars);
   TFile* fout; //< output file
   int nPars;  //< totla number of parameters
   int iStep;  //< counter for total step number
   double oldPars[NMCMCPARS]; //< array of parameters from previous step
   double oldL; //< likelihood value of previous step
   double tuneParameter;
   void setPar(int ipar,double value){oldPars[ipar]=value;}
   void setL(double value){oldL=value;}
   double varPar[NMCMCPARS]; //< stores parameter standard deviations
   void setParVar(int ipar,double value); //< sets parameter standard deviations
   void proposeStep(double* par);  //< proposes a new step from the given parameter set
   int acceptStep(double newL,double* par); 
   int acceptStepLnL(double newL,double* par);
   void savePath(const char* filename);
   void setTuneParameter(double value){tuneParameter=value;}
   TTree* pathTree;

   void test(int itry);
   TH1D* htest;
};



void markovTools::savePath(const char* filename){
  pathTree->Write();
 // fout->Close();
  pathTree->SaveAs("mcmcpath.root");
  fout->Close();
  return;
}

void markovTools::setParVar(int ipar,double value){
  varPar[ipar] = value;
  cout<<"param "<<ipar<<" variance set to: "<<value<<endl;
  return;
}

void markovTools::test(int ntry){
  nPars = 1;
  double thePar[1] = {0.};
  int itry = 0;
  int iaccept = 0;
  htest = new TH1D("htest","htest",100,-5,5);
  while (itry<ntry){
    cout<<"step: "<<itry<<endl;
    setL(1.*(thePar[0]*thePar[0]));
    cout<<"L: "<<oldL<<endl;
    proposeStep(thePar);
    cout<<"proposal: "<<thePar[0]<<endl;
    if (acceptStepLnL(1.*(thePar[0]*thePar[0]),thePar)){
      iaccept++;
      htest->Fill(thePar[0]);
    }
    cout<<"accepted: "<<iaccept<<endl; 
    itry++;
  }
  htest->Draw();
  pathTree->SaveAs("pathTree.root");
  return;
}

int markovTools::acceptStepLnL(double newL,double* par){
  double alpha = (oldL-newL);
  double rand = randy->Rndm();
  int iaccept = 0;
  //cout<<"delta L: "<<alpha<<endl;
  if (alpha>TMath::Log(rand)){
  //  pathTree->SetBranchAddress("par",par);
    for (int i=0;i<nPars;i++){
      oldPars[i]=par[i];
    }
    pathTree->Fill();
    iaccept = 1;
  } 
  else{
    for (int i=0;i<nPars;i++){
      par[i]=oldPars[i];
    }
  }
  iStep++;
  if ((iStep%100)==0) cout<<"step: "<<iStep<<endl;
  return iaccept;
}

int markovTools::acceptStep(double newL,double* par){
  double alpha = 1.;
  double lnOld;
  double lnNew;
  lnOld = TMath::Log(oldL);
  lnNew = TMath::Log(newL);
  alpha = (lnNew-lnOld);
  double rand = randy->Rndm();
  int iaccept = 0;
  if (alpha>TMath::Log(rand)){
    cout<<"accept"<<endl;
    for (int i=0;i<nPars;i++){
      oldPars[i]=par[i];
    }
    pathTree->Fill();
    iaccept = 1;
  } 
  else{
    cout<<"reject"<<endl;
    for (int i=0;i<nPars;i++){
      par[i]=oldPars[i];
    }
  }
  iStep++;
  return iaccept;
}


void markovTools::proposeStep(double* par){  
  for (int i=0;i<nPars;i++){
    //save current parameters
    oldPars[i] = par[i];
    //set new value
    par[i] = randy->Gaus(par[i],varPar[i]*tuneParameter);
   // cout<<"par: "<<oldPars[i]<<endl;
   // cout<<"parnew: "<<par[i]<<endl;
  }
  return;
}

markovTools::markovTools(int npars){
  fout = new TFile("mcmctree.root","RECREATE");
  pathTree = new TTree("MCMCpath","MCMCpath");
  nPars = npars;
  iStep = 0;
  pathTree->Branch("npars",&nPars,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",oldPars,"pars[100]/F");
}

