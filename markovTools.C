#include "TRandom2.h"
#include "TMath.h"
#include "TH1F.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#define NMCMCPARS 100

using namespace std;

TRandom2* randy = new TRandom2();

class markovTools{
   public:
   markovTools(int npars);
   int nPars;
   int iStep;
   float oldPars[NMCMCPARS];
   float oldL;
   void setPar(int ipar,float value){oldPars[ipar]=value;}
   void setL(float value){oldL=value;}
   //proposal function
   float varPar[NMCMCPARS];
   void setParVar(int ipar,float value){varPar[ipar]=value;}
   void proposeStep(float* par);
   int acceptStep(float newL,float* par);
   int acceptStepLnL(float newL,float* par);

   TTree* pathTree;

   void test(int itry);
   TH1F* htest;
};

void markovTools::test(int ntry){
  nPars = 1;
  float thePar[1] = {0.};
  int itry = 0;
  int iaccept = 0;
  htest = new TH1F("htest","htest",100,-5,5);
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

int markovTools::acceptStepLnL(float newL,float* par){
  float alpha = (oldL-newL);
  float rand = randy->Rndm();
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

int markovTools::acceptStep(float newL,float* par){
  float alpha = 1.;
  float lnOld;
  float lnNew;
  lnOld = TMath::Log(oldL);
  lnNew = TMath::Log(newL);
  alpha = (lnNew-lnOld);
  float rand = randy->Rndm();
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

void markovTools::proposeStep(float* par){  
  for (int i=0;i<nPars;i++){
    //save current parameters
    oldPars[i] = par[i];
    //set new value
    par[i] = randy->Gaus(par[i],varPar[i]);
    cout<<"par: "<<oldPars[i]<<endl;
    cout<<"parnew: "<<par[i]<<endl;
  }
  return;
}

markovTools::markovTools(int npars){
  pathTree = new TTree("MCMCpath","MCMCpath");
  nPars = npars;
  iStep = 0;
  pathTree->Branch("npars",&nPars,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",oldPars,"pars[100]/F");
}

