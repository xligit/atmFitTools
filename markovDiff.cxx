#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#include "TFile.h"

#define NMCMCPARS 200

using namespace std;

class markovDiff{
  public:

  // constructor
  markovDiff(const char* fname, int burn=0);

  // internal vars
  TFile* mcmcfile;
  TFile* outfile;
  TTree* mcmcpars;
  TTree* diffpars;
  int nburn;
  double par[NMCMCPARS];
  double par1[NMCMCPARS];
  double par2[NMCMCPARS];
  double pardiff[NMCMCPARS];
  int npars;

  // methods
  void fillDiffPars(int npairs);
  void setaddresses();
  void setouttree();

};

void markovDiff::setouttree(){
  outfile = new TFile("mcmcdiff.root","RECREATE");
  diffpars = new TTree("MCMCdiff","MCMCdiff");
  diffpars->Branch("npars",&npars,"npars/I");
  diffpars->Branch("par1",par1,"par1[200]/D");
  diffpars->Branch("par2",par1,"par2[200]/D");
  diffpars->Branch("pardiff",pardiff,"pardiff[200]/D");
  return;
}

void markovDiff::setaddresses(){
  mcmcpars->SetBranchStatus("*",0);
  mcmcpars->SetBranchStatus("npars",1);
  mcmcpars->SetBranchStatus("par",1);
  mcmcpars->SetBranchAddress("par",par);
  int nparstmp;
  mcmcpars->SetBranchAddress("npars",&nparstmp);
  mcmcpars->GetEntry(0);
  npars = nparstmp;
  cout<<"Number of mcmc pars: "<<npars<<endl;
  return; 
}

void markovDiff::fillDiffPars(int npairs){

  TRandom2* randy = new TRandom2(npairs);

  // set branch addresses of input tree
  setaddresses();

  // setup new output tree
  setouttree();

  int npts = mcmcpars->GetEntries() - nburn;

  // Loop over mcmc points and select random pairs of points after
  // some burn-in.  For these pairs, find the difference between the parameters and 
  // save this difference to the new difference tree
  for (int i=0; i<npairs; i++){
    if ((i%500)==0) cout<<i<<endl;
     cout<<i<<endl;
    // get a pair of randoms
    int ipoint = randy->Integer(npts); 
    int fpoint = randy->Integer(npts); 
    mcmcpars->GetEntry(ipoint+nburn);
    for (int ipar=0; ipar<npars; ipar++){
      par1[ipar] = par[ipar];
    }
    mcmcpars->GetEntry(fpoint+nburn);
    for (int ipar=0; ipar<npars; ipar++){
      par2[ipar] = par[ipar];
      pardiff[ipar] = par2[ipar]-par1[ipar];
    }
    diffpars->Fill(); 
  }

  diffpars->Write();
  outfile->Close();

  return;
}


markovDiff::markovDiff(const char* fname, int burn){
 
  mcmcfile = new TFile(fname);

  mcmcpars = (TTree*) mcmcfile->Get("MCMCpath");
  nburn = burn; 

  return;
}
