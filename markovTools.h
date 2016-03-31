#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#include "TFile.h"
#include "atmFitPars.h"

#define NMCMCPARS 500

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

