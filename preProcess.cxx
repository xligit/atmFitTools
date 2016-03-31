#include "preProcess.h"
#include "TObjArray.h"

////////////////////////////////////////////////////////////
//Takes in a chain and loops over all files in the chain
//For each file, a new file with a modified TTree is created
void preProcess::processAllFiles(TChain* chain){
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

void preProcess::processAllFiles(TChain *chain, TChain *banff, TChain *spline)
{
  int nfiles = chain->GetNtrees();
  TObjArray* listOfMCFiles = chain->GetListOfFiles();
  TObjArray* listOfBANFFFiles = banff->GetListOfFiles();
  TObjArray* listOfSplineFiles = spline->GetListOfFiles();
  TString tag;
  TString fmcname;
  TString fbanffname;
  TString fsplinename;
  for (int ifile=0;ifile<nfiles;ifile++){
    tag = Form("_%d_",ifile);
    fmcname = listOfMCFiles->At(ifile)->GetTitle();
    if (existBANFF) fbanffname = listOfBANFFFiles->At(ifile)->GetTitle();
    if (existSpline) fsplinename = listOfSplineFiles->At(ifile)->GetTitle();
    processFile(fmcname,fbanffname,fsplinename,tag);
  }
  return;
}

////////////////////////////////////////
//sets up pointers given a TChain
void preProcess::setTree(TChain* chin){
  TTree* trin = (TTree*)chin;
  setTree(trin);
  return;
}

///////////////////////////////////////
//sets up pointers
void preProcess::setTree(TTree* trin){
  tr = trin;
  fq = new fqReader(tr);
  vis = new visRing(fq); 
  return;
}

void preProcess::setTree(TTree* trin, TTree* banff, TTree *spline){
  tr = trin;
  trbanff = banff;
  trspline = spline;
  if (existBANFF) {
    trbanff->SetBranchAddress("fWeight", &fWeight);
    tr->AddFriend(trbanff);
    std::cout<<"setup banff tree"<<std::endl;
  }
  if (existSpline) {
    setupSplineTree(trspline);
    //tr->AddFriend(trspline);
    std::cout<<"setup spline tree "<<trspline->GetEntries()<<std::endl;
  }
  fq = new fqReader(tr);
  vis = new visRing(fq); 
  return;
}

///////////////////////////////////////////////////
// reads in a file and processes the h1 tree inside
void preProcess::processFile(const char* fname,const char* outname){

  //make output file name
  TString outputName = outDir.Data(); //name of directory
  outputName.Append(nameTag.Data()); //global name
  outputName.Append(outname); //number of this file
  outputName.Append(".root");

  //get existing h1 tree
  cout<<"opening file: "<<fname<<endl;
  TFile* fin = new TFile(fname);
  TTree* intree = (TTree*)fin->Get("h1");
  cout<<"got tree: "<<intree->GetEntries()<<endl;
  setTree(intree);

  //make new tree
  cout<<"create file: "<<outputName.Data()<<endl;
  fout = new TFile(outputName.Data(),"recreate");
  setupNewTree(); 
  std::cout<<"setup new tree"<<std::endl;
  //fill new tree
  preProcessIt();

  //clean up
  if (trout->GetEntries()>0) fout->Write();
  fin->Close();
  fout->Close();
  
  return;   
}

void preProcess::processFile(const char* fname, const char *fbanffname, const char *fsplinename, const char* outname){

  //make output file name
  TString outputName = outDir.Data(); //name of directory
  outputName.Append(nameTag.Data()); //global name
  outputName.Append(outname); //number of this file
  outputName.Append(".root");

  //get existing h1 tree
  cout<<"opening file: "<<fname<<endl;
  TFile* fin = new TFile(fname);
  TTree* intree = (TTree*)fin->Get("h1");
  cout<<"got tree: "<<intree->GetEntries()<<endl;
  TTree *banfftree=0, *splinetree=0;
  TFile *fbanffin=0, *fsplinein=0;
  if (existBANFF) {
    fbanffin = new TFile(fbanffname, "read");
    banfftree = (TTree*)fbanffin->Get("weightsTree");
    std::cout<<"opening file: "<<fbanffname<<std::endl;
    std::cout<<"got tree: "<<banfftree->GetEntries()<<std::endl;
  }
  if (existSpline) {
    fsplinein = new TFile(fsplinename, "read");
    splinetree = (TTree*)fsplinein->Get("SplinesByEvent");
    std::cout<<"opening file: "<<fsplinename<<std::endl;
    std::cout<<"got tree: "<<splinetree->GetEntries()<<std::endl;
  }

  if (splinetree->GetEntries()) setTree(intree, banfftree, splinetree);
  else setTree(intree);
  //setTree(intree);

  //make new tree
  cout<<"create file: "<<outputName.Data()<<endl;
  fout = new TFile(outputName.Data(),"recreate");
  setupNewTree(); 

  //fill new tree
  preProcessIt();

  //clean up
  if (trout->GetEntries()>0) fout->Write();
  fin->Close();
  fout->Close();
  if (fbanffin) fbanffin->Close();
  if (fsplinein) fsplinein->Close();
  return;   
}

///////////////////////////////////
//Gets a weight for an event
//Usefull for making fake data sets
float preProcess::getWeight(){
//  absmode = TMath::Abs(fq->mode);
//  float enu   = fq->pmomv[0];
  if (!existBANFF)  evtweight = 1.0;
  else if (isData)  evtweight = 1.0;
  else           evtweight = fWeight;
  //CCQE norm bin1 
//  if ((absmode==1)&&(enu<200.)){
//    evtweight = 1.5;
//  }
  //CCQE norm bin2 
//  if ((absmode==1)&&(enu>200.)&&(enu<400.)) evtweight*=1.2;
  //CCQE norm bin3 
//  if ((absmode==1)&&(enu>400.)&&(enu<800.)) evtweight*=0.9;
  //CCQE norm bin4 
//  if ((absmode==1)&&(enu>800.)) evtweight*=1.05;

  return evtweight * fq->oscwgt;
}


///////////////////////////////////////////////
//calculates the FV bin for an event
int preProcess::getBin(){

  //calculate fiducial volume variables
  //use electron hypothesis
  TVector3 vpos;
  vpos.SetXYZ(fq->fq1rpos[0][1][0],fq->fq1rpos[0][1][1],fq->fq1rpos[0][1][2]);
  TVector3 vdir;
  vdir.SetXYZ(fq->fq1rdir[0][1][0],fq->fq1rdir[0][1][1],fq->fq1rdir[0][1][2]);
  wall = calcWall(&vpos);
  towall = calcToWall(&vpos,&vdir);

  //separate into bins
  //"simple" FV Binning
  if (FVBinning==0){
    if (wall<80 && fq->fq1rmom[0][1]<1000) return 0;
    if (wall>80 && wall<200 && fq->fq1rmom[0][1]<1000) return 1;
    if (wall>200 && fq->fq1rmom[0][1]<1000) return 2;
    if (wall<80 && fq->fq1rmom[0][1]>1000) return 3;
    if (wall>80 && wall<200 && fq->fq1rmom[0][1]>1000) return 4;
    if (wall>200 && fq->fq1rmom[0][1]>1000) return 5;
  }


  //cosmic binning
//  if (FVBinning==1){
//    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
//    double Zrec = fq->fq1rpos[0][2][2];
//    double Zcut = 1410;
//    double Rcut = 1290;
//    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
//    if ((Zrec<Zcut)) return 1; //< side entering
 //   if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
//  }


  //no binning
//  if (FVBinning==2){return 0;}

  //towall binning
  if (FVBinning==3){
    if (towall<500.) return 0;
    if ((towall>-500)&&(towall<1000)) return 1;
    if (towall>=1000) return 2;
  }

  return -1;
}


/////////////////////////////////
//Simple initial cuts
int preProcess::passCuts(){
 
  //OD cut for atmospheric evetns
//  if (MCComponents==0){
//    if ((fq->nhitac)>16) return 0;
//  }

  //minimum energy cut
//  if ((fq->fqtotq[0]<80)) return 0;

  //cosmic cuts
//  if (MCComponents==1){
//    double tdecay = fq->fq1rt0[1][1]-fq->fq1rt0[0][2];
//    double ingatethresh = 1000.;
//    if ((fq->fqnse)!=2) return 0;
//    if (tdecay<ingatethresh) return 0;
//    if (fq->fq1rmom[0][1]<100.) return 0;
//    if (fq->fq1rmom[1][1]>100.) return 0;
//    if (fq->fq1rmom[1][1]<15.)  return 0;

    //towall cuts
 //   if (towall<0.) return 0; 
//  }


  //Fully Contained Cut
  if (fq->nhitac>NHITACMax) return 0;

  //Visible Energy Cut
  if (fq->fq1rmom[0][1]<EVisMin) return 0;

  //FV Basic Cuts
  if (wall<WallMin) return 0; 
  if (towall<ToWallMin) return 0;

  //Number of subevent cuts
  if (fq->fqnse>NSEMax) return 0;
  if (fq->fqnse<NSEMin) return 0;
 
  //In-gate cut
  double tdecay = fq->fq1rt0[1][1]-fq->fq1rt0[0][2];
  if (tdecay<InGateMin) return 0;


  return 1;
}

///////////////////////////////
//returns the # of subevents-1
int preProcess::getSample(){

  //atmospheric selections
  if (MCComponents==0){
    if (fq->fqnse==1) return 0;
    if (fq->fqnse==2) return 1;
    if (fq->fqnse>2)  return 2;
  }
  
  
  //cosmic selection
  if (MCComponents==1){
    return 0;
  }

  return -1;
}

int preProcess::getMode()
{
  // atm selections
  if (MCComponents==0){
    if (abs(fq->mode)==1)                          return 0; // CCQE
    else if (abs(fq->mode)>10 && abs(fq->mode)<14) return 1; // CC1pi
    else if (abs(fq->mode)==16)                    return 2; // CCCoherent
    else if (abs(fq->mode)>16 && abs(fq->mode)<30) return 3; // CCother
    else if (abs(fq->mode)>30 && abs(fq->mode)<33) return 4; // NCpi0
    else if (abs(fq->mode)>32 && abs(fq->mode)<35) return 5; // NCpipm
    else if (abs(fq->mode)==36)                    return 6; // NCCoherent
    else if (abs(fq->mode)>36 && abs(fq->mode)<53) return 7; // NCother (including NCel)
    else if (abs(fq->mode)==2)                     return 8; // MEC
  }  
  
  //cosmic selection
  if (MCComponents==1){
    return 0;
  }

  return -1;
}

//////////////////////////////////////
//get code for MC true component type
int preProcess::getComponent(){

  //useful for cuts
  absmode = TMath::Abs(fq->mode);
  int absnu   = TMath::Abs(fq->ipnu[0]);
 
  //visible event selection
  //charged current
  if (MCComponents==0){
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
  }


  //cosmic selectoin
  if (MCComponents==1){
    return 0;
  }

}

//loop over all events and sort into bins, samples and components
void preProcess::preProcessIt(){
  int nev = tr->GetEntries();
  for (int i=0;i<nev;i++){
    //get info for event
    tr->GetEntry(i);
    if (!isData) {if (nFilesSpline>1) trspline->GetEntry(i);}
    if ((i%1000)==0) cout<<i<<endl;
    nbin=getBin();
    if (!passCuts()) continue;
    vis->fillVisVar(); //get visible ring information
    fillAttributes();
    ncomponent=getComponent();
    nsample=getSample();
    nmode=getMode();
    evtweight=getWeight();
    trout->Fill();
  }

  return;
}

///////////////////////////////////////
//returns the index of the best 2R fit
int preProcess::getBest2RFitID(){
  
  int nfits = fq->fqnmrfit;

  double ngLnLBest = 1000000.;
  int bestindex = 0;

  for (int ifit=0;ifit<nfits;ifit++){
    int fitID = fq->fqmrifit[ifit]; //< fit fit ID code
    if ((fitID-320000000)<0) continue; //< we want best 2R fits
    if (fq->fqmrnll[ifit]<ngLnLBest) bestindex = ifit;
  }
  return bestindex;
}

////////////////////////////////////////
//fills fiTQun attribute array
void preProcess::fillAttributes(){
  attribute[0] = fq->fq1rnll[0][2]-fq->fq1rnll[0][1];
//  if (fq->fqnse>1) attribute[1] = fq->fq1rnll[1][2]-fq->fq1rnll[1][1];
//  else{
//    attribute[1] = 0.;
//  }
  int ibest = getBest2RFitID();
  double best1Rnglnl = fmin(fq->fq1rnll[0][1],fq->fq1rnll[0][2]);
  attribute[1] = best1Rnglnl-fq->fqmrnll[ibest];

//  attribute[3] = -fq->fqmrnll[ibest];
//  attribute[2] = fq->fq1rmom[0][1];
//  attribute[3] = fq->fq1rmom[0][2];
//  attribute[4] = fq->fq1rmom[1][1];
//  if (fq->fqnse>1){
//    attribute[1] = fq->fq1rnll[1][2]-fq->fq1rnll[1][1];
//    attribute[4] = fq->fq1rmom[1][1];
//  }
//  else{
//    attribute[1] = 0.;
//    attribute[4] = 0.;
//  }
  return;
}

void preProcess::setupNewTree(){
  tr->SetBranchStatus("*",0);
  tr->SetBranchStatus("fq*",1);
  tr->SetBranchStatus("*v",1);
  tr->SetBranchStatus("ipnu",1);
  //tr->SetBranchStatus("cluster*",1);
  tr->SetBranchStatus("mode",1);
  tr->SetBranchStatus("nhitac",1);
  tr->SetBranchStatus("oscwgt");
  trout = tr->CloneTree(0); //clone but don't copy data
  trout->CopyAddresses(tr); //set addresses
  
  //add new branches
  trout->Branch("attribute",attribute,"attribute[1000]/F");
  trout->Branch("ncomponent",&ncomponent,"ncomponent/I");
  trout->Branch("nsample",&nsample,"nsample/I");
  trout->Branch("nmode",&nmode,"nmode/I");
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
  if (!isData) {
    trout->Branch("byEv_maqe_ccqe_gr", "TGraph", &byEv_maqe_ccqe_gr, 1280000, 0);
    trout->Branch("byEv_pfo_ccqe_gr", "TGraph", &byEv_pfo_ccqe_gr, 1280000, 0);
    trout->Branch("byEv_ebo_ccqe_gr", "TGraph", &byEv_ebo_ccqe_gr, 1280000, 0);
    trout->Branch("byEv_ca5_cc1pi_gr", "TGraph", &byEv_ca5_cc1pi_gr, 1280000, 0);
    trout->Branch("byEv_ca5_ncpiz_gr", "TGraph", &byEv_ca5_ncpiz_gr, 1280000, 0);
    trout->Branch("byEv_ca5_ncpipm_gr", "TGraph", &byEv_ca5_ncpipm_gr, 1280000, 0);
    trout->Branch("byEv_manff_cc1pi_gr", "TGraph", &byEv_manff_cc1pi_gr, 1280000, 0);
    trout->Branch("byEv_manff_ncpiz_gr", "TGraph", &byEv_manff_ncpiz_gr, 1280000, 0);
    trout->Branch("byEv_manff_ncpipm_gr", "TGraph", &byEv_manff_ncpipm_gr, 1280000, 0);
    trout->Branch("byEv_bgscl_cc1pi_gr", "TGraph", &byEv_bgscl_cc1pi_gr, 1280000, 0);
    trout->Branch("byEv_bgscl_ncpiz_gr", "TGraph", &byEv_bgscl_ncpiz_gr, 1280000, 0);
    trout->Branch("byEv_bgscl_ncpipm_gr", "TGraph", &byEv_bgscl_ncpipm_gr, 1280000, 0);
    trout->Branch("byEv_dismpishp_ccoth_gr", "TGraph", &byEv_dismpishp_ccoth_gr, 1280000, 0);
    //trout->Branch("byEv_sccvec_ccqe_gr", "TGraph",  &byEv_sccvec_ccqe_gr, 1280000, 0);
    //trout->Branch("byEv_sccvec_ncoth_gr", "TGraph",  &byEv_sccvec_ncoth_gr, 1280000, 0);
    //trout->Branch("byEv_sccaxl_ccqe_gr", "TGraph",  &byEv_sccaxl_ccqe_gr, 1280000, 0);
    //trout->Branch("byEv_sccaxl_ncoth_gr", "TGraph",  &byEv_sccaxl_ncoth_gr, 1280000, 0);
    trout->Branch("byEv_rpa_ccqe_gr", "TGraph",  &byEv_rpa_ccqe_gr, 1280000, 0);
  }
  return;
}

void preProcess::setupSplineTree(TTree *h)
{
  h->SetBranchAddress("byEv_maqe_ccqe_gr", &byEv_maqe_ccqe_gr, &byEv_maqe_ccqe_br);
  h->SetBranchAddress("byEv_pfo_ccqe_gr", &byEv_pfo_ccqe_gr, &byEv_pfo_ccqe_br);
  h->SetBranchAddress("byEv_ebo_ccqe_gr", &byEv_ebo_ccqe_gr, &byEv_ebo_ccqe_br);
  h->SetBranchAddress("byEv_ca5_cc1pi_gr", &byEv_ca5_cc1pi_gr, &byEv_ca5_cc1pi_br);
  h->SetBranchAddress("byEv_ca5_ncpiz_gr", &byEv_ca5_ncpiz_gr, &byEv_ca5_ncpiz_br);
  h->SetBranchAddress("byEv_ca5_ncpipm_gr", &byEv_ca5_ncpipm_gr, &byEv_ca5_ncpipm_br);
  h->SetBranchAddress("byEv_manff_cc1pi_gr", &byEv_manff_cc1pi_gr, &byEv_manff_cc1pi_br);
  h->SetBranchAddress("byEv_manff_ncpiz_gr", &byEv_manff_ncpiz_gr, &byEv_manff_ncpiz_br);
  h->SetBranchAddress("byEv_manff_ncpipm_gr", &byEv_manff_ncpipm_gr, &byEv_manff_ncpipm_br);
  h->SetBranchAddress("byEv_bgscl_cc1pi_gr", &byEv_bgscl_cc1pi_gr, &byEv_bgscl_cc1pi_br);
  h->SetBranchAddress("byEv_bgscl_ncpiz_gr", &byEv_bgscl_ncpiz_gr, &byEv_bgscl_ncpiz_br);
  h->SetBranchAddress("byEv_bgscl_ncpipm_gr", &byEv_bgscl_ncpipm_gr, &byEv_bgscl_ncpipm_br);
  h->SetBranchAddress("byEv_dismpishp_ccoth_gr", &byEv_dismpishp_ccoth_gr, &byEv_dismpishp_ccoth_br);
  h->SetBranchAddress("byEv_sccvec_ccqe_gr", &byEv_sccvec_ccqe_gr, &byEv_sccvec_ccqe_br);
  h->SetBranchAddress("byEv_sccvec_ncoth_gr", &byEv_sccvec_ncoth_gr, &byEv_sccvec_ncoth_br);
  h->SetBranchAddress("byEv_sccaxl_ccqe_gr", &byEv_sccaxl_ccqe_gr, &byEv_sccaxl_ccqe_br);
  h->SetBranchAddress("byEv_sccaxl_ncoth_gr", &byEv_sccaxl_ncoth_gr, &byEv_sccaxl_ncoth_br);
  h->SetBranchAddress("byEv_rpa_ccqe_gr", &byEv_rpa_ccqe_gr, &byEv_rpa_ccqe_br);
}

/////////////////////////////
//empty constructor
preProcess::preProcess(){
  nFiles=0;
  
  byEv_maqe_ccqe_gr = 0;
  byEv_pfo_ccqe_gr = 0;
  byEv_ebo_ccqe_gr = 0;
  byEv_ca5_cc1pi_gr = 0;
  byEv_ca5_ncpiz_gr = 0;
  byEv_ca5_ncpipm_gr = 0;
  byEv_manff_cc1pi_gr = 0;
  byEv_manff_ncpiz_gr = 0;
  byEv_manff_ncpipm_gr = 0;
  byEv_bgscl_cc1pi_gr = 0;
  byEv_bgscl_ncpiz_gr = 0;
  byEv_bgscl_ncpipm_gr = 0;
  byEv_dismpishp_ccoth_gr = 0;
  byEv_sccvec_ccqe_gr = 0;
  byEv_sccvec_ncoth_gr = 0;
  byEv_sccaxl_ccqe_gr = 0;
  byEv_sccaxl_ncoth_gr = 0;
  byEv_rpa_ccqe_gr = 0;

}

/////////////////////////////
//construct from TTree
preProcess::preProcess(TTree* trin,const char* name){
  tr=trin;
  fq = new fqReader(tr);
  vis = new visRing(fq);
  setupNewTree();
  nameTag=name;
  return;
}

/////////////////////////////
//construct from TChain
preProcess::preProcess(TChain* chin,const char* name){
  tr=(TTree*)chin;
  nameTag=name;
  fq = new fqReader(tr);
  vis = new visRing(fq);
  setupNewTree();
  return;
}

preProcess::preProcess(TChain *mc, TChain *banff, TChain *spline, const std::string name)
{
  tr = (TTree*)mc;
  trbanff = (TTree*)banff;
  if (trbanff) trbanff->SetBranchAddress("fWeight", &fWeight);
  trspline = (TTree*)spline;
  setupSplineTree(trspline);
  nameTag = name.c_str();
  fq = new fqReader(tr);
  vis = new visRing(fq);

  setupNewTree();
}

//////////////////////////////////////////
//read in parameters and run preprocessing!
void preProcess::runPreProcessing(){
  
  //read in parameters!
  cout<<"preProcess: Reading parameters from file "<<parFileName.Data()<<endl;
  sharedPars* runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
  nameTag = runpars->globalRootName;
  cout<<"nametag: "<<nameTag.Data()<<endl;
  FVBinning = runpars->preProcessFVBinning; //< flag for FV binning type in getBin()
  MCComponents = runpars->preProcessMCComponents; //< flag for MC component definitions in getComponent()
  NHITACMax = runpars->PreProcFCCut;
  EVisMin = runpars->PreProcEVisCut;
  WallMin = runpars->PreProcWallMinCut;
  ToWallMin = runpars->PreProcToWallMinCut;
  NSEMax = runpars->PreProcNseMax0;
  NSEMin = runpars->PreProcNseMin;
  InGateMin = runpars->PreProcInGateCut; 

  //create data and mc chains
  chmc = new TChain("h1");
  chdat = new TChain("h1");
  chbanff = new TChain("weightsTree");
  chspline = new TChain("SplinesByEvent");
  chmc->Add(runpars->preProcessFilesMC.Data());
  chdat->Add(runpars->preProcessFilesData.Data());
  existBANFF = false; existSpline = false;
  nFilesBANFF = chbanff->Add(runpars->preProcessFilesBANFF.Data());
  if (nFilesBANFF) existBANFF = true;
  nFilesSpline = chspline->Add(runpars->preProcessFilesSpline.Data());
  if (nFilesSpline) existSpline = true;
  std::cout<<"# banff files = "<<nFilesBANFF<<"; # spline files = "<<nFilesSpline<<std::endl;
  if (chmc->GetEntries()<1){
    cout<<"preProcess ERROR: no events in MC chain"<<endl;
    return;
  }
  if (chdat->GetEntries()<1){
    cout<<"preProcess ERROR: no events in Data chain"<<endl;
    return;
  }

  //process the files
  outDir = runpars->preProcessOutDir.Data();
  nameTag.Append("_ppmc");
  isData = false;
  if (existBANFF || existSpline) processAllFiles(chmc, chbanff, chspline);
  else processAllFiles(chmc);
  nameTag = runpars->globalRootName.Data();
  nameTag.Append("_ppdata");
  isData = true;
  processAllFiles(chdat); 

  cout<<"preProcess: Complete!"<<endl;

  //////////////////////////
  return;
}


