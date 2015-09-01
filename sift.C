//#include "/nfs/hepusers/users/amissert/stdROOTinc.h"
#include <iostream>
#include "visRing.C"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"




using namespace std;



class sift{
  public:
  sift(TChain* chin);
  sift(TTree* trin);
  TTree* tr;
  TTree* trout;
  fqReader* fq;
  visRing*  vis;
  void setupNewTree();
  void siftIt(const char* filename);
  int ncomponent;
  int nsample;
  int nbin;
  //selections
  int absmode;
  int getComponent();
  int getSample();
  int getBin();
};

int sift::getBin(){
  return 0;
}

int sift::getSample(){
  if (fq->fqnse==1) return 1;
  if (fq->fqnse==2) return 2;
  if (fq->fqnse>2)  return 3;
  return 0;
}

int sift::getComponent(){
  absmode = TMath::Abs(fq->mode);
  if ((absmode)<30){
     if ((vis->nve==1)&&(vis->nvpip+vis->nvpi0+vis->nvk)==0) return 1; //CC1e
     if ((vis->nvmu==1)&&(vis->nvpip+vis->nvpi0+vis->nvk)==0) return 2; //CC1mu
     if (TMath::Abs(fq->ipnu[0])==14) return 4; //CCmuOther
     if (TMath::Abs(fq->ipnu[0])==12) return 3; //CCeOther
  }
  else{
     if ((vis->nvpi0==1)&&(vis->nvpip==0)) return 5;
     return 6;
  }
  return 7;
}

void sift::siftIt(const char* fname){
  int nev = tr->GetEntries();
  for (int i=0;i<nev;i++){
    //get info for event
    if ((i%1000)==0) cout<<i<<endl;
    tr->GetEntry(i);
    vis->fillVisVar();
    //sort event
    ncomponent=getComponent();
    nsample=getSample();
    nbin=getBin();
    trout->Fill();
  }

  TString name = fname;
  name.Append(".root");
  trout->SaveAs(name.Data());
  return;
}

void sift::setupNewTree(){
  trout = tr->CloneTree(0); //clone but don't copy data
  trout->CopyAddresses(tr); //set addresses
  
  //add new branches
  trout->Branch("ncomponent",&ncomponent,"ncomponent/I");
  trout->Branch("nsample",&nsample,"nsample/I");
  trout->Branch("nbin",&nbin,"nbin/I");
  trout->Branch("nvis",&vis->nvis,"nvis/I");
  trout->Branch("nvmu",&vis->nvmu,"nvmu/I");
  trout->Branch("nve",&vis->nve,"nve/I");
  trout->Branch("nvgam",&vis->nvgam,"nvgam/I");
  trout->Branch("nvpip",&vis->nvpip,"nvpip/I");
  trout->Branch("nvpi0",&vis->nvpi0,"nvpi0/I");
  trout->Branch("nvoth",&vis->nvoth,"nvoth/I");
  return;
}

sift::sift(TTree* trin){
  tr=trin;
  fq = new fqReader(tr);
  vis = new visRing(fq);
  setupNewTree();
  return;
}

sift::sift(TChain* chin){
  tr=(TTree*)chin;
  fq = new fqReader(tr);
  vis = new visRing(fq);
  setupNewTree();
  return;
}



