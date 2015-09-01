#include "/nfs/hepusers/users/amissert/stdROOTinc.h"
#include <iostream>
#include "visRing.C"
#include "fqReaderCosmic.C"

using namespace std;

TRandom2 *randy = new TRandom2();

class tweaker{
  public:
  tweaker(TTree* trin);
  tweaker(TChain* chin);
  TTree* tr;
  TTree* trout;
  fqReaderCosmic* fq;
  void Init();

  double Alpha;
  double Beta;
  double Gamma;
  double meanDecayMom;

  double tweakDecayMom1();
  double tweakDecayMom2(double mean);
  double tDecayMom1;
  double tDecayMom2;
  void makeTweakTree(const char* nametag);
};

void tweaker::makeTweakTree(const char* nametag){
  TString fname = nametag;
  fname.Append(".root");
//  TTree* trout = tr->CloneTree(0);
//  trout->CopyAddresses(tr);
  int nev=tr->GetEntries();
  for (int i=0;i<nev;i++){
    tr->GetEntry(i);
    tDecayMom1=tweakDecayMom1();
    tDecayMom2=tweakDecayMom2(meanDecayMom);
    cout<<"tDecayMom2: "<<tDecayMom2<<endl;
    trout->Fill();
  }
  trout->SaveAs(fname.Data());
  return;
}

double tweaker::tweakDecayMom2(double mean){
  double delta = randy->Gaus(0.,Gamma);  
  return Alpha*(fq->fq1rmom[1][1]+delta-mean)+mean+Beta;
}

double tweaker::tweakDecayMom1(){
  double delta = randy->Gaus(Beta,Alpha);
  return fq->fq1rmom[1][1]+delta;
}

void tweaker::Init(){
  fq = new fqReaderCosmic(tr);
  tr->SetBranchStatus("*",0);
  tr->SetBranchStatus("fq*",1);
  trout = tr->CloneTree(0);
  trout->CopyAddresses(tr);
  trout->Branch("tDecayMom1",&tDecayMom1,"tDecayMom1/D");
  trout->Branch("tDecayMom2",&tDecayMom2,"tDecayMom2/D");
  cout<<"tweaker got made."<<endl;
}

tweaker::tweaker(TChain* chin){
  tr = (TTree*)chin;
  Init();
  return;
}


tweaker::tweaker(TTree* trin){
  tr = trin;
  Init();
  return;
}

