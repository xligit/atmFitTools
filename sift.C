#ifndef SIFT_C
#define SIFT_C

#include "sift.h"
#include "TObjArray.h"

void sift::processAllFiles(TChain* chain){
  int nfiles = chain->GetNtrees();
  TObjArray* listOfFiles = chain->GetListOfFiles();
  TString tag;
  TString fname;
  for (int ifile=0;ifile<nfiles;ifile++){
    tag = Form("_%d_",ifile);
    fname = listOfFiles->At(ifile)->GetTitle();
    processFile(fname,tag);
  }
  return;
}

TString sift::getFileRootName(){
  return nameTag;
}

void sift::setTree(TChain* chin){
  TTree* trin = (TTree*)chin;
  setTree(trin);
  return;
}

void sift::setTree(TTree* trin){
/*  if (vis!=NULL){
     cout<<"delete visible ring counter"<<endl;
     delete vis;
     vis=NULL;
  }
  if (fq!=NULL){
     cout<<"deleting fiTQun reader"<<endl;
     delete fq;
     fq=NULL;
  }
  */
  tr = trin;
  fq = new fqReader(tr);
  vis = new visRing(fq); 
  return;
}

void sift::processFile(const char* fname,const char* outname){

  //make output file name
  TString outputName = nameTag.Data();
  outputName.Append(outname);
  outputName.Append(".root");
//  outputName.Append(fname);
  //

  //get tree
  cout<<"opening file: "<<fname<<endl;
  TFile* fin = new TFile(fname);
  TTree* intree = (TTree*)fin->Get("h1");
  cout<<"got tree: "<<intree->GetEntries()<<endl;
  setTree(intree);

  //make new tree
  cout<<"create file: "<<outputName.Data()<<endl;
  fout = new TFile(outputName.Data(),"recreate");
  setupNewTree(); 

  //fill new tree
  siftIt();
//
  //delete new tree
  //trout->Delete();
  fout->Write();
  fin->Close();
  fout->Close();
  
  return;   
}

void sift::addFile(const char*filename){
  fileNames[nFiles]=filename;
  nFiles++;
}

float sift::getWeight(){
  int absmode = TMath::Abs(fq->mode);
  float enu   = fq->pmomv[0];
  evtweight = 1.0;
  //CCQE norm bin1 
  if ((absmode==1)&&(enu<200.)){
//    evtweight = 1.5;
  }
  //CCQE norm bin2 
//  if ((absmode==1)&&(enu>200.)&&(enu<400.)) evtweight*=1.2;
  //CCQE norm bin3 
//  if ((absmode==1)&&(enu>400.)&&(enu<800.)) evtweight*=0.9;
  //CCQE norm bin4 
//  if ((absmode==1)&&(enu>800.)) evtweight*=1.05;

  return evtweight;
}

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
  //charged current
  if ((absmode>0)&&(absmode<30)){
    if ((vis->nve==1)&&((vis->nvp==0)&&(vis->nvmu==0)&&(vis->nvpi0==0)&&(vis->nvpip==0))) return 0; //CC single e
    if ((vis->nvmu==1)&&((vis->nvp==0)&&(vis->nvpi0==0)&&(vis->nvpip==0))) return 1; //CC single mu
    if ((vis->nve==1)&&(vis->nvmu==0)) return 2; //CC e + other
    if ((vis->nve==0)&&(vis->nvmu==1)) return 3; //CC mu + other
    return 4; //CC Other
  }
  else{
    if ((vis->nvpi0>0)) return 5;  //single pi0
    return 6;  //NC other
  }
//  return 7;


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
void sift::siftIt(){
  int nev = tr->GetEntries();
  for (int i=0;i<nev;i++){
    //get info for event
    if ((i%1000)==0) cout<<i<<endl;
    tr->GetEntry(i);
    if (!passCuts()) continue;
    vis->fillVisVar(); //get visible ring information
    ncomponent=getComponent();
    nsample=getSample();
    nbin=getBin();
    evtweight=getWeight();
    trout->Fill();
  }
 // TString name = fname;
 // name.Append(".root");
//  trout->SaveAs(name.Data());
 // trout->Write();
 // fout->Write();
 // fout->Close();
 // trout->Delete();
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
//  TString filename=nameTag.Data();
//  filename.Append("_siftOutput.root");
//  fout = new TFile(filename.Data(),"recreate");
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
  trout->Branch("evtweight",&evtweight,"evtweight/F");
  return;
}

sift::sift(const char* name){
  nameTag=name;
  nFiles=0;
}

sift::sift(TTree* trin,const char* name){
  tr=trin;
  fq = new fqReader(tr);
  vis = new visRing(fq);
  setupNewTree();
  nameTag=name;
  return;
}

sift::sift(TChain* chin,const char* name){
  tr=(TTree*)chin;
  nameTag=name;
  fq = new fqReader(tr);
  vis = new visRing(fq);
  setupNewTree();
  return;
}


#endif
