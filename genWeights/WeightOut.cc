
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <TString.h>
#include <TArrayF.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>

using namespace std;

int main(int argc, char* argv[]){
  
  TString ifname=argv[1];
  Float_t *fWeight;
  TArrayF *weights= new TArrayF();
  Int_t nweights,nwgt;
  
  TFile ofile(ifname+"_weights.root","RECREATE");
  TChain trIn("weightstree");
  trIn.Add(ifname+"_wgt.root");
  trIn.SetBranchAddress("weights", &weights);
  trIn.SetBranchAddress("nweights",&nweights);
  trIn.GetEntry(0);
  nwgt=nweights;
  fWeight = new Float_t[nwgt];
  
  TTree fOutputTree("weightstree", "Tree storing TArray of weights and twkdials per event");
  fOutputTree.Branch("fWeight", &fWeight[0], Form("fWeight[%d]/F",nwgt));
  
  cout << " # of entries: " << trIn.GetEntries() << endl;
  for (int i=0; i<trIn.GetEntries(); i++) {
    trIn.GetEntry(i);
    for (int j=0; j<nwgt; j++) {
      fWeight[j]=weights->fArray[j];
    }
    fOutputTree.Fill();
  }
  
  fOutputTree.Write();
  ofile.Close();
  
  delete weights;
  delete[] fWeight;
  
}
