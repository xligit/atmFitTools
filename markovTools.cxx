#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#include "TFile.h"
#include "atmFitPars.cxx"

#define NMCMCPARS 500

using namespace std;

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
#ifdef T2K
TRandom3 *randy = new TRandom3();
#else
TRandom2* randy = new TRandom2();
#endif
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
   void setL(double value){oldL=value; cout<<"Lset: "<<value<<endl;}
   void setParVar(int ipar,double value); //< sets parameter standard deviations

   /////////////////////////
   //MCMC Functions
   void proposeStep(double* par);  //< proposes a new step from the given parameter set
   void proposeStep(); //< propose a step from the parameters in atmFitPars
   int  acceptStepLnL(double newL); //< decide if step is accepted given new LnL
   int acceptStep(double newL,double* par); 
   int acceptStepLnL(double newL,double* par);

   /////////////////////////
   //I/O
   void savePath();
   void setTuneParameter(double value);

   /////////////////////////
   //debugging
   void test(int itry);
   void testAtm(int itry);

   TH1D* htest;
   TH2D* htest2D;


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
  cout<<"parameter "<<ipar<<" variance set to: "<<value<<endl;
  return;
}

void markovTools::testAtm(int ntry){

  //setup atmFitPars
  atmFitPars* fitpars = new atmFitPars(1,1,1,1);
  fitpars->setParameter(0,0.);
  fitpars->setParameter(1,0.);
  fitpars->setSysParUnc(0,1.);
  fitpars->setSysParUnc(1,1.);
  nPars = fitpars->nTotPars;
  atmPars = fitpars;

  tuneParameter = 1.0;

  int itry = 0;
  int iaccept = 0;
  double LnL;
  htest2D = new TH2D("htest2","htest2",50,-5,5,50,-5,5);
  while (itry<ntry){
    cout<<"step: "<<itry<<endl;
    double par1 = atmPars->getParameter(0);
    double par2 = atmPars->getParameter(1);
    LnL = par1*par1 + par2*par2;
    setL(LnL);
    cout<<"L: "<<oldL<<endl;
    proposeStep();
    par1 = atmPars->getParameter(0);
    par2 = atmPars->getParameter(1);
    LnL = par1*par1 + par2*par2;
    if (acceptStepLnL(LnL)){
      iaccept++;
      htest2D->Fill(atmPars->getParameter(0),atmPars->getParameter(1));
    }
    cout<<"accepted: "<<iaccept<<endl; 
    itry++;
  }
  htest2D->Draw("colz");
  //pathTree->SaveAs("pathTree.root");
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

/////////////////////////////////////////////
//decide if new parameters shoudl be accepted
int markovTools::acceptStepLnL(double newL){


  double alpha = (oldL-newL); //< get difference in LnL
 // cout<<"mcmc: Likelihood: "<<oldL;;
//  cout<<" -> "<<newL;
//  cout<<" diff: "<<alpha<<endl;
  double rand = randy->Rndm(); //< throw a random number
 // cout<<"Log unit random: "<<TMath::Log(rand)<<endl;
  int iaccept = 0; //< acceptance flag

  /////////////////////////////
  //chekc if we should accept
  if (alpha>TMath::Log(rand)){
    //accepted! new pars are now old
    for (int i=0;i<nPars;i++){
      oldPars[i]=atmPars->getParameter(i);
      oldL = newL;
    }
#ifdef T2K
    atmPars->acceptStep();
#endif
    pathTree->Fill();
    iaccept = 1;
  } 
  else{
    for (int i=0;i<nPars;i++){
      //rejected, reset to old parameters
      atmPars->setParameter(i,oldPars[i]);
    }
  }

//  if (iaccept) cout<<"accepted!"<<endl;
//  else{
//   cout<<"not accepted"<<endl;
//  }
  iStep++; //< increment global step  count
  if ((iStep%100)==0) cout<<"step: "<<iStep<<endl;

  /////////////////////////////
  return iaccept;
}



int markovTools::acceptStepLnL(double newL,double* par){
  double alpha = (oldL-newL);
  cout<<"oldL: "<<oldL<<endl;
  cout<<"newL: "<<newL<<endl;
  cout<<"Ldiff: "<<alpha<<endl;
  double rand = randy->Rndm();
  cout<<"rand: "<<rand<<endl;
  int iaccept = 0;
  //cout<<"delta L: "<<alpha<<endl;
  if (alpha>TMath::Log(rand)){
  //  pathTree->SetBranchAddress("par",par);
    for (int i=0;i<nPars;i++){
      oldPars[i]=par[i];
    }
#ifdef T2K
    atmPars->acceptStep();
#endif
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
  if (iaccept) cout<<" accepted! "<<endl;
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
#ifdef T2K
  atmPars->proposeStep();
#endif
  for (int i=0;i<nPars;i++){
#ifndef T2K
    oldPars[i] = par[i];
    if (fixPar[i]!=1){
      par[i] = randy->Gaus(par[i],varPar[i]*tuneParameter);
      cout<<"Proposed parameter "<<i<<" as :"<<par[i]<<endl;
    }
#else
    if (fixPar[i]!=1) par[i] = atmPars->getPropParameter(i);
#endif
  }
  return;
}


/////////////////////////////////////////////////
//takes a poiter to atmFitPars and suggests new parameters
//
void markovTools::proposeStep(){  
#ifndef T2K
  for (int i=0;i<nPars;i++){
    oldPars[i] = atmPars->getParameter(i);
    if (atmPars->fixPar[i]!=1){
//      cout<<"par "<<i<<": "<<oldPars[i];
      double random = randy->Gaus(oldPars[i],atmPars->parUnc[i]*tuneParameter);
      //cout<<"random: "<<random<<endl;
      atmPars->setParameter(i,random);
 //     cout<<"->  "<<atmPars->getParameter(i)<<endl;;
    }
  }
#else
  atmPars->proposeStep();
  for (int i = 0; i < nPars; ++i) {
    oldPars[i] = atmPars->getParameter(i);
  }
#endif
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
  pathTree->Branch("par",oldPars,"par[500]/D");

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
  pathTree->Branch("par",oldPars,"par[500]/D");
}

////////////////////////////////////////////
//construct from atmFitPars 
markovTools::markovTools(atmFitPars* fitpars){
  fout = new TFile("mcmctree.root","RECREATE"); //< set output file name
  pathTree = new TTree("MCMCpath","MCMCpath"); //< initialize new tree for steps
  atmPars = fitpars;
  Init(fitpars->nTotPars); 
  return;
}

void markovTools::setTuneParameter(double value)
{
  tuneParameter=value;
#ifdef T2K
  atmPars->setStepSize(tuneParameter);
#endif
}
