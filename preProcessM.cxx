#ifndef PREPROCESS_C
#define PREPROCESS_C

#include "preProcess.h"
#include "TObjArray.h"

/////////////////////////////////////////////////////////////////
// Setup the FV bin histogram for getFVBin()
void preProcess::setFVBinHisto(){
  hFVBins = new TH2FV("hfvbins",1);
}


/////////////////////////////////////////////////////////////////
//Create a TGraph that is used to make weights to see effect of
//correcting for cosmic flux mis-modeling
void  preProcess::setWeightHistogram(const char* file, const char * name){

  //read in histogram
  TFile* hfile = new TFile(file);
  hWeight = (TH1D*)hfile->Get(name);
 
  //convert to graph
  const int nn = hWeight->GetNbinsX();
  double xx[nn];
  double yy[nn];
  for (int i=0;i<nn;i++){
    xx[i] = hWeight->GetBinCenter(i+1);
    yy[i] = hWeight->GetBinContent(i+1);
  }
  gWeight = new TGraph(nn,xx,yy);

  gWeight->Draw("alc");

  useWeights = 1;

  return;
}

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

void preProcess::processAllFiles(TChain *chain, TChain *spline)
{
  int nfiles = chain->GetNtrees();
  TObjArray* listOfMCFiles = chain->GetListOfFiles();
  TObjArray* listOfSplineFiles = spline->GetListOfFiles();
  TString tag;
  TString fmcname;
  TString fsplinename;
  for (int ifile=0;ifile<nfiles;ifile++){
    tag = Form("_%d_",ifile);
    fmcname = listOfMCFiles->At(ifile)->GetTitle();
    if (existSpline) fsplinename = listOfSplineFiles->At(ifile)->GetTitle();
    processFile(fmcname,fsplinename,tag);
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
  fq = new fqEvent(tr);
  vis = new visRing(fq); 
  return;
}

void preProcess::setTree(TTree* trin, TTree *spline){
  tr = trin;
  trspline = spline;
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
  int neventsnew = preProcessIt();
  cout<<"  filled it with "<<neventsnew<<" events!" <<endl;

  //clean up
  //if (trout->GetEntries()>0) {
  fout->cd();
  trout->Write();
  //}
  fin->Close();
  fout->Close();
  
  return;   
}

void preProcess::processFile(const char* fname, const char *fsplinename, const char* outname){

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
  TTree *splinetree=0;
  TFile *fsplinein=0;
  if (existSpline) {
    std::cout<<"processing spline"<<std::endl;
    fsplinein = new TFile(fsplinename, "read");
    splinetree = (TTree*)fsplinein->Get("SplinesByEvent");
    std::cout<<"opening file: "<<fsplinename<<std::endl;
    std::cout<<"got tree: "<<splinetree->GetEntries()<<std::endl;
  }

  if (splinetree->GetEntries()) setTree(intree, splinetree);
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
  if (fsplinein) fsplinein->Close();
  return;   
}

///////////////////////////////////
//Gets a weight for an event
//Usefull for making fake data sets
float preProcess::getBANFFWeight(){
//  absmode = TMath::Abs(fq->mode);
//  float enu   = fq->pmomv[0];
  float weight = 1;
  if (isData)  weight = 1.0;
  else  {
    weight = fq->fWeight;
    weight *= fq->wgtflx[3] * fq->wgtosc1[3];
  }
  return weight;
}



void preProcess::setParFileName(const char* fname){
  
  // set file name
  parFileName=fname;
  
  //read in parameters!
  cout<<"preProcess: Reading parameters from file "<<parFileName.Data()<<endl;
  sharedPars* runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
  nameTag = runpars->globalRootName;
  cout<<"nametag: "<<nameTag.Data()<<endl;
  FVBinning = runpars->preProcessFVBinning; //< flag for FV binning type in getBin()
  setFVBinHisto();
  MCComponents = runpars->preProcessMCComponents; //< flag for MC component definitions in getComponent()
  MCSamples = runpars->preProcessMCSamples; //< flag for MC sample definitions in getSample()
  NHITACMax = runpars->preProcFCCut;
  EVisMin = runpars->preProcEVisCut;
  WallMin = runpars->preProcWallMinCut;
  ToWallMin = runpars->preProcToWallMinCut;
  NSEMax = runpars->preProcNseMax0;
  NSEMin = runpars->preProcNseMin;
  InGateMin = runpars->preProcInGateCut; 
  flgAddMoreVars = runpars->preProcAddMoreVars;

  // list of attributes to use
  nAttributes = runpars->nAttributes;
  attributeList[0] = runpars->fQAttName0;
  attributeList[1] = runpars->fQAttName1;
  attributeList[2] = runpars->fQAttName2;
  attributeList[3] = runpars->fQAttName3;
  attributeList[4] = runpars->fQAttName4;
  attributeList[5] = runpars->fQAttName5;
  attributeList[6] = runpars->fQAttName6;
  attributeList[7] = runpars->fQAttName7;


float preProcess::getFixedWeight()
{
  float weight = 1;
  if (isData)  weight = 1.0;
  else  {
    weight = fq->rfgWeight;
    weight *= fq->wgtflx[3] * fq->wgtosc1[3];
  }
  return weight;
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
  // calculate additional fv variables as well
  for (int isubev=0; isubev<fq->fqnse; isubev++){
    vpos.SetXYZ(fq->fq1rpos[isubev][2][0],fq->fq1rpos[isubev][2][1],fq->fq1rpos[isubev][2][2]);
    vdir.SetXYZ(fq->fq1rdir[isubev][2][0],fq->fq1rdir[isubev][2][1],fq->fq1rdir[isubev][2][2]);
    fq1rwall[isubev][2] = calcWall2(&vpos);
    fq1rtowall[isubev][2] = calcToWall(&vpos,&vdir);
    vpos.SetXYZ(fq->fq1rpos[isubev][1][0],fq->fq1rpos[isubev][1][1],fq->fq1rpos[isubev][1][2]);
    vdir.SetXYZ(fq->fq1rdir[isubev][1][0],fq->fq1rdir[isubev][1][1],fq->fq1rdir[isubev][1][2]);
    fq1rwall[isubev][1] = calcWall2(&vpos);
    fq1rtowall[isubev][1] = calcToWall(&vpos,&vdir);
  }
  //true towall
  for (int ipart=0; ipart<fq->npar; ipart++){
    vpos.SetXYZ(fq->posv[0],fq->posv[1],fq->posv[2]);
    vdir.SetXYZ(fq->dirv[ipart][0],fq->dirv[ipart][1],fq->dirv[ipart][2]);
    towallv[ipart]=calcToWall(&vpos,&vdir);
    wallv2=calcWall2(&vpos);
  }
  // even more FV related variables
  if (flgAddMoreVars>0){
     // add 1r perimiter and mincone	  
     for (int isubev=0; isubev<fq->fqnse; isubev++){
       vpos.SetXYZ(fq->fq1rpos[isubev][2][0],fq->fq1rpos[isubev][2][1],fq->fq1rpos[isubev][2][2]);
       vdir.SetXYZ(fq->fq1rdir[isubev][2][0],fq->fq1rdir[isubev][2][1],fq->fq1rdir[isubev][2][2]);
       fq1rperim[isubev][2] = calcPerimeter(&vpos,&vdir);
       fq1rmincone[isubev][2] = calcMinCone(&vpos,&vdir);
       vpos.SetXYZ(fq->fq1rpos[isubev][1][0],fq->fq1rpos[isubev][1][1],fq->fq1rpos[isubev][1][2]);
       vdir.SetXYZ(fq->fq1rdir[isubev][1][0],fq->fq1rdir[isubev][1][1],fq->fq1rdir[isubev][1][2]);
       fq1rperim[isubev][1] = calcPerimeter(&vpos,&vdir);
       fq1rmincone[isubev][1] = calcMinCone(&vpos,&vdir);
    }
    //true perimeter and mincone
    if (flgAddMoreVars>1){
      for (int ipart=0; ipart<fq->npar; ipart++){
        vpos.SetXYZ(fq->posv[0],fq->posv[1],fq->posv[2]);
        vdir.SetXYZ(fq->dirv[ipart][0],fq->dirv[ipart][1],fq->dirv[ipart][2]);
        perimv[ipart]=calcPerimeter(&vpos,&vdir);
        minconev[ipart]=calcMinCone(&vpos,&vdir);
      }
    }
  }

  ////////////////////////////////
  // separate into bins by FVBinning parameter

  //////////////////////////////////////////
  //"simple" FV Binning for atm
  if (FVBinning==0){
    if (wall<80 && fq->fq1rmom[0][1]<1000) return 0;
    if (wall>80 && wall<200 && fq->fq1rmom[0][1]<1000) return 1;
    if (wall>200 && fq->fq1rmom[0][1]<1000) return 2;
    if (wall<80 && fq->fq1rmom[0][1]>1000) return 3;
    if (wall>80 && wall<200 && fq->fq1rmom[0][1]>1000) return 4;
    if (wall>200 && fq->fq1rmom[0][1]>1000) return 5;
  }

  ////////////////////////////////////////
  // cosmic binning by entering surface
  if (FVBinning==1){
    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
    double Zrec = fq->fq1rpos[0][2][2];
    double Zcut = 1410;
    double Rcut = 1290;
    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
    if ((Zrec<Zcut)) return 1; //< side entering
    if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
  }

  //////////////////////////////////
  // all in one bin
  if (FVBinning==2){return 0;}
  
  //////////////////////////////////
  // towall binning for cosmics
  if (FVBinning==3){
    if (towall<500.) return 0;
    if ((towall>-500)&&(towall<1000)) return 1;
    if (towall>=1000) return 2;
  }

  //////////////////////////////////////
  // binning using wall/towallhistogram
  if (FVBinning==4){
      int fvbin = hFVBins->FindBin(towall,wall)-1;
      return fvbin;
   }
//    if (fvbin<0){
//      cout<<"Bad FV value:"<<endl;
//      cout<<"  wall:   "<<wall<<endl;
//      cout<<"  towall: "<<towall<<endl;
//  //  }
      //return fvbin;
 //   if (fq->fq1rmom[0][1]<230.) return fvbin;
  //  if (fq->fq1rmom[0][1]>800.) return fvbin+12;
   // return fvbin+6;
  

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

  //std::cout<<fq->nhitac<<" "<<fq->fq1rmom[0][1]<<" "<<wall<<" "<<towall<<" "<<fq->fqnse<<"\n"
  //	   <<NHITACMax<<" "<<EVisMin<<" "<<WallMin<<" "<<ToWallMin<<" "<<NSEMax<<" "<<NSEMin<<"\n"<<std::endl;

  //Fully Contained Cut
  if ((int)fq->nhitac>NHITACMax) return 0;

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

  //hybrid pi0 selection
  if (MCSamples==2){
//    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
//    double Zrec = fq->fq1rpos[0][2][2];
//    double Zcut = 1410;
//    double Rcut = 1290;
//    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
//    if ((Zrec<Zcut)) return 1; //< side entering
//    if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
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
  //std::cout<<"nev = "<<nev<<std::endl;
  for (int i=0;i<nev;i++){
    //get info for event
    tr->GetEntry(i);
    if (!isData) {
      if (existSpline) trspline->GetEntry(i);
    }
    if ((i%1000)==0) cout<<i<<endl;
    nbin=getBin();
    if (!passCuts()) continue;
    vis->fillVisVar(); //get visible ring information
    fillAttributes();
    ncomponent=getComponent();
    nsample=getSample();
    nmode=getMode();
    evtweight = getBANFFWeight();
    rfgweight = getFixedWeight();
    trout->Fill();
    //if(isData) std::cout<<i<<" filling data "<<std::endl;
  }
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
  tr->SetBranchStatus("oscwgt",1);
  tr->SetBranchStatus("wgtosc1",1);
  tr->SetBranchStatus("wgtosc2",1);
  tr->SetBranchStatus("wgtflx",1);
  tr->SetBranchStatus("fWeight",1);
  tr->SetBranchStatus("rfgWeight",1);
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
  trout->Branch("rfgweight",&rfgweight,"rfgweight/F");

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
    //trout->Branch("byEv_rpa_ccqe_gr", "TGraph",  &byEv_rpa_ccqe_gr, 1280000, 0);
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
  //h->SetBranchAddress("byEv_sccvec_ccqe_gr", &byEv_sccvec_ccqe_gr, &byEv_sccvec_ccqe_br);
  //h->SetBranchAddress("byEv_sccvec_ncoth_gr", &byEv_sccvec_ncoth_gr, &byEv_sccvec_ncoth_br);
  //h->SetBranchAddress("byEv_sccaxl_ccqe_gr", &byEv_sccaxl_ccqe_gr, &byEv_sccaxl_ccqe_br);
  //h->SetBranchAddress("byEv_sccaxl_ncoth_gr", &byEv_sccaxl_ncoth_gr, &byEv_sccaxl_ncoth_br);
  //h->SetBranchAddress("byEv_rpa_ccqe_gr", &byEv_rpa_ccqe_gr, &byEv_rpa_ccqe_br);
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

preProcess::preProcess(TChain *mc, TChain *spline, const std::string name)
{
  tr = (TTree*)mc;
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
  chspline = new TChain("SplinesByEvent");
  chmc->Add(runpars->preProcessFilesMC.Data());
  chdat->Add(runpars->preProcessFilesData.Data());
  existSpline = false;
  nFilesSpline = chspline->Add(runpars->preProcessFilesSpline.Data());
  if (nFilesSpline) existSpline = true;
  std::cout<<"# spline files = "<<nFilesSpline<<std::endl;
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
  if (existSpline) processAllFiles(chmc, chspline);
  else processAllFiles(chmc);
  nameTag = runpars->globalRootName.Data();
  nameTag.Append("_ppdata");
  isData = true;
  std::cout<<"---------- process data -----------"<<std::endl;
  processAllFiles(chdat); 

  cout<<"preProcess: Complete!"<<endl;

  //////////////////////////
  return;
}


