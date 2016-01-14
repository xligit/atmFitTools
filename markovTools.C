#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#include "TFile.h"
#include "atmFitPars.C"

#define NMCMCPARS 100

using namespace std;

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
TRandom2* randy = new TRandom2();
#endif

//class to manage a Markov Chain Monte Carlo
class markovTools{
   public:
   
   ///////////////////
   //constructors
   markovTools(int npars);
   markovTools(atmFitPars* atmpars);
   void Init(int pars);

   ///////////////////////////
   //variables
   TFile* fout; //< output file
   int nPars;  //< totla number of parameters
   int iStep;  //< counter for total step number
   double oldPars[NMCMCPARS]; //< array of parameters from previous step
   int fixPar[NMCMCPARS]; //< array of fix flags for each parameter
   double oldL; //< likelihood value of previous step
   double tuneParameter; //< tunes the size of MCMC steps
   double varPar[NMCMCPARS]; //< stores parameter standard deviations
   TTree* pathTree;
   atmFitPars* atmPars;

   /////////////////////////
   //setters
   void setFixPar(int ipar, int value){fixPar[ipar]=value;}
   void setPar(int ipar,double value){oldPars[ipar]=value;}
   void setL(double value){oldL=value;}
   void setParVar(int ipar,double value); //< sets parameter standard deviations

   /////////////////////////
   //MCMC Functions
   void proposeStep(double* par);  //< proposes a new step from the given parameter set
   int acceptStep(double newL,double* par); 
   int acceptStepLnL(double newL,double* par);

   /////////////////////////
   //I/O
   void savePath();
   void setTuneParameter(double value){tuneParameter=value;}

   /////////////////////////
   //debugging
   void test(int itry);
   TH1D* htest;

};



void markovTools::savePath(){
  pathTree->Write();
 // fout->Close();
//  if (filename){
//    pathTree->SaveAs(filename);
//  }
//  else{
//    pathTree->SaveAs("mcmcpath.root");
//  }
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



/////////////////////////////////////////////////
//takes a poiter to a parameter array and returns
//a new array of proposed parameters
void markovTools::proposeStep(double* par){  
  for (int i=0;i<nPars;i++){
    oldPars[i] = par[i];
    if (fixPar[i]!=1) par[i] = randy->Gaus(par[i],varPar[i]*tuneParameter);
  }
  return;
}

/////////////////////////////////////////////
//initialization
void markovTools::Init(int npars){

  // set total number of parameters
  nPars = npars;

  //current step counter
  iStep=0;
  
  //branch setup
  pathTree->Branch("npars",&nPars,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",oldPars,"pars[100]/F");

  //done
  return;
}

/////////////////////////////////////////////
//constructor
markovTools::markovTools(int npars){
  fout = new TFile("mcmctree.root","RECREATE"); //< set output file name
  pathTree = new TTree("MCMCpath","MCMCpath"); //< initialize new tree for steps
  nPars = npars; //< set total # of parameters
  iStep = 0;  //< counter of current step
  
  //////////////////////////////////////////
  //branch setup
  pathTree->Branch("npars",&nPars,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",oldPars,"pars[100]/F");
}

////////////////////////////////////////////
//construct from atmFitPars 
markovTools::markovTools(atmFitPars* fitpars){
  Init(fitpars->nTotPars);

  return;
}



