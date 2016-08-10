#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"

#include "shared.h"

#include "atmFitPars.cxx"
#include "histoManager.cxx"
#include "fqProcessedEvent.cxx"
#include "splineFactory.cxx"
#include "getSystWeight.cxx"
#include "mcmcReader.cxx"
#include "THStack.h"
#include "sharedPars.cxx"
#include "TH2FV.h"


/////////////////////////////////////////////////////////
// Class to manage and apply parameters from the MCMC fit
class markovParManager{
  public:
  
  /////////////////////////////////////////////////////////////
  // construct from a parameter file
  markovParManager(const char* parfile);

  void setMCMCFile(const char* filename);

  atmFitPars* fitpars; //< manages the fit parameters
  sharedPars* runpars; //< runtime parameters
  TFile* mcmcfile;
  TTree* mcmcpath;
  // leafs
  double mcmcpar[500];
  int nmcmcpar;

  void setParsToPoint(int ipoint); //< sets all parameters to a point in mcmc cloud

  void modifyAttributes(fqProcessedEvent* fqevent, int mcmcpt); //< sets parameters and modifies attributes 

  void 

};


///////////////////////////////////////////////////////////////////
// modify the attribues for a given event
void markovParManager::modifyAttributes(fqProcessedEvent* fqevent, int mcmcpt){
  int nattributes = fitPars->nAttributes;
  
  int parindex = fitPars->getParIndex(  
  double biaspar = fitPars->getPar 

   for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
    int evtbin = fqevent->nbin;
    int evtcomp =  fqevent->ncomponent;
    int evtsamp = fqevent->nsample;
    double mean = hManager->hMCMean[evtsamp][evtbin][evtcomp][iatt];
    double scaling = fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0];
    double bias    = fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1];
    mcreader->attribute[iatt] *= scaling; //multiply by smear parameter
    mcreader->attribute[iatt] += (mean - (scaling*mean)); //< correct mean shift from scaling   
    mcreader->attribute[iatt] += bias; //add bias parameter
  }

  //////////
  return;


}

////////////////////////////////////////////////////////////////////
// set parameters in internal parameter arrays to a given mcmc point
void markovParManager::setParsToPoint(int ipoint){
  for (int ipar=0; ipar<nmcmcpar; ipar++){
    fitPars->setParameter(ipar,mcmcpar[ipar]);
  }
  return;
}


////////////////////////////////////////////////////////
// reads in the mcmc path tree from a file
void markovParManager::setMCMCFile(const char* filename){
   mcmcfile = new TFile(filename);

   mcmcpath = (TTree*)mcmcfile->Get("MCMCpath");

   // set addresses
   partree->SetBranchAddress("par",mcmcpar);
   partree->SetBranchAddress("npars",&nmcmcpar);
   return;
}


markovParManager::markovParManager(const char* parfile){

  // read in parameters
  runPars = new sharedPars(parfile);
  runPars->readParsFromFile();
 
  fitpars = new atmFitPars(parfile);


}



