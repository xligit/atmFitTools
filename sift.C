#ifndef SIFT_C
#define SIFT_C

#include "sift.h"


int sift::getBin(){
  //calculate fiducial volume variables
  //use electron hypothesis
  TVector3 vpos;
  vpos.SetXYZ(fq->fq1rpos[0][1][0],fq->fq1rpos[0][1][1],fq->fq1rpos[0][1][2]);
  TVector3 vdir;
  vdir.SetXYZ(fq->fq1rdir[0][1][0],fq->fq1rdir[0][1][1],fq->fq1rdir[0][1][2]);
  wall = calcWall(&vpos);
  towall = calcToWall(&vpos,&vdir);
  if ((wall<200.)&&(wall>50.)) return 1;
  if (wall<50.) return 2;
  return 0;
}

int sift::passCuts(){
  if ((fq->nhitac)>16) return 0;
  if ((fq->fqtotq[0]<80)) return 0;
  return 1;
}

int sift::getSample(){
  if (fq->fqnse==1) return 0;
  if (fq->fqnse==2) return 1;
  if (fq->fqnse>2)  return 2;
  return 3;
}

int sift::getComponent(){
  absmode = TMath::Abs(fq->mode);
  int absnu   = TMath::Abs(fq->ipnu[0]);
  if ((vis->nve==1)&&((vis->nvk==0)&&(vis->nvmu==0)&&(vis->nvpi0==0)&&(vis->nvpip==0))) return 0; //single e
  if ((vis->nvmu==1)&&((vis->nvk==0)&&(vis->nvpi0==0)&&(vis->nvpip==0))) return 1; //single mu
  if ((vis->nve==1)&&(vis->nvmu==0)&&((vis->nvpi0>0)||(vis->nvpip>0)||(vis->nvk>0))) return 2; //e + other
  if ((vis->nve==0)&&(vis->nvmu==1)&&((vis->nvpi0>0)||(vis->nvpip>0)||(vis->nvk>0))) return 3; //mu + other
  if ((vis->nve==0)&&(vis->nvmu==0)&&(vis->nvpi0>0)&&(vis->nvpip==0)) return 4;  //single pi0
  if ((vis->nve==0)&&(vis->nvmu==0)&&(vis->nvpi0==0)&&(vis->nvpip>0)) return 5;  //single pip
 // if ((vis->nvpi0>0)||(vis->nvpip>0)||(vis->nvk>0)) return 6;  //OTHER
  return 6;


//  if ((absmode)<30){
//     if ((absmode==1)&&(absnu==14)) return 1;
//     if ((absmode==1)&&(absnu==12)) return 0;
  //   if ((vis->nvpip==0)&&(vis->nvpip+vis->nvpi0+vis->nvk)==0) return 1; //CC1e
 //    if ((vis->nvmu==1)&&(vis->nvpip+vis->nvpi0+vis->nvk)==0) return 2; //CC1mu
//     if (TMath::Abs((int)fq->ipnu[0])==14) return 3; //CCmuOther
//     if (TMath::Abs((int)fq->ipnu[0])==12) return 2; //CCeOther
//  }
//  else{
//     if ((vis->nvpi0==1)&&(vis->nvpip==0)) return 4;
//     return 5;
//   }
//  return 6;
}

//loop over all events and sort into bins, samples and components
void sift::siftIt(const char* fname){
  int nev = tr->GetEntries();
  for (int i=0;i<nev;i++){
    //get info for event
    if ((i%1000)==0) cout<<i<<endl;
    tr->GetEntry(i);
    if (!passCuts()) continue;
    vis->fillVisVar(); //get visible ring information
    //sort event
    ncomponent=getComponent();
  //  if (ncomponent==3) fq->fq1rnll[0][1]-=100.;
 //   if (ncomponent==2) fq->fq1rnll[0][2]+=(-1*25.);
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
  tr->SetBranchStatus("*",0);
  tr->SetBranchStatus("fq*",1);
  tr->SetBranchStatus("*v",1);
  tr->SetBranchStatus("ipnu",1);
  //tr->SetBranchStatus("cluster*",1);
  tr->SetBranchStatus("mode",1);
  tr->SetBranchStatus("nhitac",1);
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
  trout->Branch("nvp",&vis->nvp,"nvp/I");
  trout->Branch("nvk",&vis->nvk,"nvk/I");
  trout->Branch("nvoth",&vis->nvoth,"nvoth/I");
  trout->Branch("vispid",vis->vispid,"vispid[100]/I");
  trout->Branch("wall",&wall,"wall/F");
  trout->Branch("towall",&towall,"towall/F");
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


#endif
