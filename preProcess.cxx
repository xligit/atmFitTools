#ifndef PREPROCESS_C
#define PREPROCESS_C

#include "preProcess.h"

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

  // find how many files are in this chain
  int nfiles = chain->GetNtrees();
  
  // get a list of all of the files
  TObjArray* listOfFiles = chain->GetListOfFiles();

  // loop over files and process each one
  TString tag;
  TString fname;
  for (int ifile=0;ifile<nfiles;ifile++){
    tag = Form("_%d",ifile);
    fname = listOfFiles->At(ifile)->GetTitle();
    processFile(fname,tag);
  }

  //
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
  fq = new fqEvent(tr, ntupleType.Data());
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
  cout<<"  opening file: "<<fname<<endl;
  TFile* fin = new TFile(fname);
  TTree* intree = (TTree*)fin->Get("h1");
  cout<<"  got tree: "<<intree->GetEntries()<<endl;
  setTree(intree); //< set pointers to current tree

  //make new tree
  cout<<"  create file: "<<outputName.Data()<<endl;
  fout = new TFile(outputName.Data(),"recreate");
  setupNewTree(); 

  //fill new tree
  int neventsnew = preProcessIt();
  cout<<"  filled it with "<<neventsnew<<" events!" <<endl;


  //clean up
  if (neventsnew>0) fout->Write();
  fin->Close();
  fout->Close();
  
  return;   
}

///////////////////////////////////////////////////////////////////
// makes a couple test bumps in the specified directory
void preProcess::makeTestFiles(const char* outdir, int testtype, int nmc, int ndata, int randseed){

  // random numbers
  TRandom2* randy = new TRandom2(randseed);

  // file name
  TString fnamebase = "testbump_";

  // make a simple gaussian bump
  if (testtype==0){
   
   // fill fake mc files
   outDir = outdir;
   TString outputName = outDir.Data(); //name of directory
   outputName.Append(fnamebase.Data()); 
   outputName.Append("_ppmc.root"); 

   // make mc file
   fout = new TFile(outputName.Data(),"recreate");

   //make new mc tree and set branches
   trout = new TTree("h1","h1");
   trout->Branch("attribute",attribute,"attribute[1000]/F");
   trout->Branch("fqrcpar",&fqrcpar,"fqrcpar/F");
   trout->Branch("ncomponent",&ncomponent,"ncomponent/I");
   trout->Branch("nsample",&nsample,"nsample/I");
   trout->Branch("nbin",&nbin,"nbin/I");

   // fill with fake data
   for (int i=0; i<nmc; i++){
      attribute[0] = randy->Gaus();
      nsample = 0;
      nbin    = 0;
      ncomponent = 0;
      trout->Fill();
   }

   // write output
   fout->Write();
   fout->Close();
  
   // now fill fake data file
   outputName = outDir.Data(); //name of directory
   outputName.Append(fnamebase.Data()); 
   outputName.Append("_ppdata.root"); 

   // make data file
   fout = new TFile(outputName.Data(),"recreate");

   //make new tree and set branches
   trout = new TTree("h1","h1");
   trout->Branch("attribute",attribute,"attribute[1000]/F");
   trout->Branch("fqrcpar",&fqrcpar,"fqrcpar/F");
   trout->Branch("ncomponent",&ncomponent,"ncomponent/I");
   trout->Branch("nsample",&nsample,"nsample/I");
   trout->Branch("nbin",&nbin,"nbin/I");

   // fill with fake data
   for (int i=0; i<ndata; i++){
      attribute[0] = randy->Gaus();
      nsample = 0;
      nbin    = 0;
      ncomponent = 0;
      trout->Fill();
   }

   // write output
   fout->Write();
   fout->Close();
  }


  return;
}

///////////////////////////////////
//Gets a weight for an event
//Usefull for making fake data sets
float preProcess::getWeight(){
  evtweight = 1.0;
  if (useWeights){
    evtweight = gWeight->Eval(fq->fq1rmom[0][2],0,"s");
  }

// if using skimmed tree from Xiaoyue, calculate event weight
// based on these variables in the ntuples
#ifdef USE_Xl_WEIGHTS


  evtweight *= fq->wgtosc1[3]; // for Xiaoyue MC
  evtweight *= fq->wgtflx[3]; // for Xiaoyue MC

#endif
#ifdef USE_ST_WEIGHTS

  //////////////////////////////////////////////
  // For Shimpei MC
  double w_maxsolact =0.35;
  evtweight *= ((1. - w_maxsolact)*fq->flxh11[0]+w_maxsolact*fq->flxh11[2])/fq->flxh06[1];
  evtweight *= fq->oscwgt; // for Shimpei MC

#endif

  // fake normalization bump
  if (fakeNormFlg){
    if (fq->pmomv[0]<1000.){
      evtweight*=1.15;
    }
  }

  return evtweight;
}



void preProcess::setParFileName(const char* fname){
  
  // set file name
  parFileName=fname;
   
  //read in parameters!
  cout<<"preProcess: Reading parameters from file "<<parFileName.Data()<<endl;
  sharedPars* runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
  nameTag = runpars->globalRootName;
  FVBinning = runpars->preProcessFVBinning; //< flag for FV binning type in getBin()
  if (FVBinning==4) setFVBinHisto();
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
  flgUseSpikeMask = runpars->preProcMaskFlg;
  if (flgUseSpikeMask>0){
     TString fmaskname = runpars->preProcMaskFile.Data();
     TFile* maskfile = new TFile(fmaskname.Data());
     cout<<"preProcess: Getting spike mask from file: "<<fmaskname.Data()<<endl;     
     hmask = (TH1D*)maskfile->Get("hmask");
  }

  // type of ntuple structure in input files
  // used to avoid printed errors in "SetBranchAddress"
  ntupleType = runpars->ntupleType;
  cout<<"preProcess: ntuple type: "<<ntupleType.Data()<<endl;

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


  //create data and mc chains
  chmc = new TChain("h1");
  chdat = new TChain("h1");
  cout<<"preProcess: adding MC files: "<<runpars->preProcessFilesMC.Data()<<endl;
  cout<<"preProcess: adding Data files: "<<runpars->preProcessFilesData.Data()<<endl;
  chmc->Add(runpars->preProcessFilesMC.Data());
  chdat->Add(runpars->preProcessFilesData.Data());
  if (chmc->GetEntries()<1){
    cout<<"preProcess WARNING: no events in MC chain"<<endl;
  }
  if (chdat->GetEntries()<1){
    cout<<"preProcess WARNING: no events in Data chain"<<endl;
  }

  outDir = runpars->preProcessOutDir.Data();
  


  return;
  /*
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
  */

}



///////////////////////////////////////////////
//calculates the FV bin for an event
int preProcess::getBin(){

  ////////////////////////////////////////////////////////////
  //calculate fiducial volume variables
  //use electron hypothesis
  TVector3 vpos;
  vpos.SetXYZ(fq->fq1rpos[0][2][0],fq->fq1rpos[0][2][1],fq->fq1rpos[0][2][2]);
  TVector3 vdir;
  vdir.SetXYZ(fq->fq1rdir[0][2][0],fq->fq1rdir[0][2][1],fq->fq1rdir[0][2][2]);
  wall = calcWall2(&vpos);
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
    if ((wall<200.)&&(wall>80.)) return 1;
    if (wall<80.) return 0;
    return 2;
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

  

  return -1;
}


/////////////////////////////////
//Simple initial cuts
int preProcess::passCuts(){

  /////////////////////
  //tmp cuts
 // if (towallv[0]<80.) return 0;

  /////////////////////
  //Fully Contained Cut
  if ((int)fq->nhitac>NHITACMax) return 0;

  //////////////////////
  //Visible Energy Cut
  if (fq->fq1rmom[0][1]<EVisMin) return 0;

  ////////////////
  //FV Basic Cuts
  if (wall<WallMin) return 0; 
  if (towall<ToWallMin) return 0;  

  /////////////////////////
  //Number of subevent cuts
  if (fq->fqnse>NSEMax) return 0;
  if (fq->fqnse<NSEMin) return 0;
 
  /////////////
  //In-gate cut
  if (InGateMin>0){
    double tdecay = fq->fq1rt0[1][1]-fq->fq1rt0[0][2];
    if (tdecay<InGateMin) return 0;
  }

  /////////////////////////////////////////////
  // optional masking cut for hybrid pi0 spikes
  if (flgUseSpikeMask>0){
     if (!passMask(hmask,fq1rwall[0][2])) return 0;
  }

  ////////////////////
  //all cuts passed!!
  return 1;
}

///////////////////////////////
//returns the # of subevents-1
int preProcess::getSample(){

  //atmospheric selections
  if (MCSamples==0){
    if (fq->fqnse==1) return 0;
    if (fq->fqnse==2) return 1;
    if (fq->fqnse>2)  return 2;
  }

  //atmospheric selection with energy
  if (MCSamples==3){
    double evis = fq->fq1rmom[0][1];
    if (evis<1000.){
      if (fq->fqnse==1) return 0;
      if (fq->fqnse==2) return 1;
      if (fq->fqnse>2)  return 2;
    }
    else{
      if (fq->fqnse==1) return 3;
      if (fq->fqnse==2) return 4;
      if (fq->fqnse>2)  return 5;
    }
  }

  
  //cosmic selection
  if (MCSamples==1){
//    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
//    double Zrec = fq->fq1rpos[0][2][2];
//    double Zcut = 1410;
//    double Rcut = 1290;
//    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
//    if ((Zrec<Zcut)) return 1; //< side entering
//    if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
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

  cout<<"preProcess:getSample():  Warning sample code"<<MCSamples<<" not defined!"<<endl;
  return -1;
}

//////////////////////////////////////
//get code for MC true component type
int preProcess::getComponent(){

  ////////////////////////////
  // useful for cuts
  absmode = TMath::Abs(fq->mode);
 // int absnu   = TMath::Abs(fq->ipnu[0]);
 
  /////////////////////////////////////////
  // visible + NEUT event selection for atm
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

  /////////////////////////////////////////
  // visible only component selections
  if (MCComponents==2){

    ///////////////////////////////////////////
    // 0 -> Single Shower
    // 1 -> Single MIP
    // 2 -> Shower + other (not single pi0)
    // 3 -> MIP + other
    // 4 -> Single Pi0
    // 5 -> Other (should be zero events)
    //////////////////////////////////////////

    double showerthresh = 15.;
    double nonshowerthresh = 45.;

    // no visible (decay e)
    if (vis->nvis==0){
      return 0;
    }
   
    // weak 1R (decay e)
    if (vis->nvis==1 && vis->visbrightness[0]<nonshowerthresh) return 0;

    // single pi0
    if (vis->nvis<=2){
       // all visible rings are gammas and there is single pi0
       if ((vis->nvis-vis->nvgam==0) && (vis->nvpi0==1)) return 4;
    }

    // select single ring events 
    if (vis->nvis==1){
      if (vis->nve==1) return 0;
      if (vis->nvgam==1) return 0;
      if (vis->nvmu==1) return 1;
      if (vis->nvpip==1) return 1;
      if (vis->nvp==1)   return 1;
    }

    // select weak MR events
    if (vis->nvis>1 &&
       ( (vis->vismrbrightness<showerthresh && vis->vismrtype2==1) ||
         (vis->vismrbrightness<nonshowerthresh && vis->vismrtype2==0) )){
      if (vis->vismrtype1==0) return 1; //< count as non-showering ring
      if (vis->vismrtype1==1) return 0; //< count as showering ring
    }

    // select MR events
    if (vis->nvis>1){
      // showering most visible ring
      if (vis->vismrtype1 == 1) return 2;
      // non-showering most visible ring;
      if (vis->vismrtype1 == 0) return 3;
    }

    // should be nothing?
    else{
      return 5;
    }
    

  };
 
  //////////////////////////////////////////
  // cosmic selectoin
  if (MCComponents==1){
    return 0;
  }
  

  //////////////////////////////////////////
  // hybrid pi0 selectoin
  if (MCComponents==3){
    return 0;
  }

  cout<<"preProcess: ERROR MC component not defined!";
  return -1;
}

//loop over all events and sort into bins, samples and components
int preProcess::preProcessIt(){
  int nev = tr->GetEntries();
  int naccepted = 0;
  for (int i=0;i<nev;i++){
    //get info for event
    if ((i%1000)==0) cout<<"event:  "<<i<<endl;

    tr->GetEntry(i);
    //calc FV bin and fill FV variables
    nbin=getBin();
    if (nbin<0.) continue;
    //apply cuts
    if (!passCuts()) continue;
    naccepted++;
    // hybrid pi0s don't have the right banks for VR counting
    if (MCComponents!=3) vis->fillVisVar(); //get visible ring information
    fillAttributes(fq);
    ncomponent=getComponent();
    if (fakeShiftFlg){
      if (ncomponent==0){
//        attributeMap["fqelike"] = attributeMap["fqelike"] + 50.;
        attribute[0] = attribute[0] +50.;
      }
    }
    nsample=getSample();
    evtweight=getWeight();
    trout->Fill();
  }

  return naccepted;
}

///////////////////////////////////////
//returns the index of the best 2R fit
int preProcess::getBest2RFitID(){
  
  int nfits = fq->fqnmrfit;

  double ngLnLBest = 10000000.;
  int bestindex = 0;

  for (int ifit=0;ifit<nfits;ifit++){
    int fitID = TMath::Abs(fq->fqmrifit[ifit]); //< fit fit ID code
    if ((TMath::Abs(fitID-20000000))>100) continue; //< we want best 2R fits
    if (fq->fqmrnll[ifit]<ngLnLBest){
      bestindex = ifit;
      ngLnLBest=fq->fqmrnll[ifit];
    }
  }
  best2RID = fq->fqmrifit[bestindex];
  return bestindex;
}

////////////////////////////////////////
//fills fiTQun attribute array
void preProcess::fillAttributes(fqEvent* fqevent){
;
  // Fill the cmap that matches attribute names to values
  fillAttributeMap(fqevent);

  // fill the attribute[] array with the values you want to use for this analysis
  for (int i=0; i<nAttributes; i++){
    double value = attributeMap[attributeList[i].Data()];
    attribute[i] = value;
  }

  return;
}


////////////////////////////////////////
//fills cmap of possible fitqun attributes
void preProcess::fillAttributeMap(fqEvent* fqevent){

  // PID e vs. mu ratio of first subevent
  attributeMap["fqelike"] = fqevent->fq1rnll[0][2]-fqevent->fq1rnll[0][1];

  // Ring Counting (RC) parameter
  int ibest = getBest2RFitID();
  double best1Rnglnl = fmin(fqevent->fq1rnll[0][1],fqevent->fq1rnll[0][2]);
  fqrcpar = best1Rnglnl-fqevent->fqmrnll[ibest];
  attributeMap["fqrcpar"] = fqrcpar;

  // Reconstructed distance from wall
  attributeMap["fqwall"] = wall;

  // Reconstructed momentum (muon)
  attributeMap["fq1rmumom"] = fqevent->fq1rmom[0][2];

  // Reconstructed momentum (electron)
  attributeMap["fq1remom"] = fqevent->fq1rmom[0][1];

  // Reconstructed momentum (best 1R)
  if ((fqevent->fq1rnll[0][2]-fqevent->fq1rnll[0][1])>0.){
    attributeMap["fq1rmom"] = fqevent->fq1rmom[0][1];
  }
  else{
    attributeMap["fq1rmom"] = fqevent->fq1rmom[0][2];
  }

  // pi0 likelihood
  attributeMap["fqpi0like"] = fqevent->fq1rnll[0][1] - fqevent->fqpi0nll[0];

  // pi0 mass
  attributeMap["fqpi0mass"] = fqevent->fqpi0mass[0];

  // pi0 photon angle
  attributeMap["fqpi0photangle"] = fqevent->fqpi0photangle[0];

  // pi0 wall min
  TVector3 vpos;
  vpos.SetXYZ(fqevent->fqpi0pos[0][0],fqevent->fqpi0pos[0][1],fqevent->fqpi0pos[0][2]);
  TVector3 vdir1;
  vdir1.SetXYZ(fqevent->fqpi0dir1[0][0],fqevent->fqpi0dir1[0][1],fqevent->fqpi0dir1[0][2]);
  TVector3 vdir2;
  vdir2.SetXYZ(fqevent->fqpi0dir2[0][0],fqevent->fqpi0dir2[0][1],fqevent->fqpi0dir2[0][2]);
  double pi0wall= calcWall2(&vpos);
  double pi0towall1 = calcToWall(&vpos,&vdir1);
  double pi0towall2 = calcToWall(&vpos,&vdir2);
  attributeMap["fqpi0wall"] = pi0wall;

  // pi0 towall min
  attributeMap["fqpi0towallmin"] = fmin(pi0towall1,pi0towall2);

  // pi0 total momentum
  attributeMap["fqpi0momtot"] = fqevent->fqpi0momtot[0];

  return;
}

void preProcess::setupNewTree(){

  // select just a few branches
  tr->SetBranchStatus("*",0);
  tr->SetBranchStatus("fq*",1);
  tr->SetBranchStatus("*v",1);
  tr->SetBranchStatus("ipnu",1);
  tr->SetBranchStatus("mode",1);
  tr->SetBranchStatus("nhitac",1);
  tr->SetBranchStatus("nring",1);
  tr->SetBranchStatus("*scnd*",1);
  tr->SetBranchStatus("nscndprt",1);
  tr->SetBranchStatus("iprnttrk",1);
  tr->SetBranchStatus("iprntprt",1);
  tr->SetBranchStatus("iorgprt",1);
  tr->SetBranchStatus("iprntidx",1);
  tr->SetBranchStatus("nchilds",1);
  tr->SetBranchStatus("ichildidx",1);
  tr->SetBranchStatus("iflgscnd",1);
#ifdef USE_XL_WEIGHTS
  tr->SetBranchStatus("wgt*",1);
#endif
#ifdef USE_ST_WEIGHTS
  tr->SetBranchStatus("oscwgt",1);
  tr->SetBranchStatus("flxh11",1);
  tr->SetBranchStatus("flxh06",1);
#endif
  tr->SetBranchStatus("*vc",1);

  // make new output tree that is clone of old tree,
  // with a few extra branchs
  trout = tr->CloneTree(0); //clone but don't copy data
  trout->CopyAddresses(tr); //set addresses
  //add new branches
  trout->Branch("attribute",attribute,"attribute[1000]/F");
  trout->Branch("fqrcpar",&fqrcpar,"fqrcpar/F");
  trout->Branch("ncomponent",&ncomponent,"ncomponent/I");
  trout->Branch("nsample",&nsample,"nsample/I");
  trout->Branch("nbin",&nbin,"nbin/I");
  // visible ring counting
  trout->Branch("nvis",&vis->nvis,"nvis/I");
  trout->Branch("nvmu",&vis->nvmu,"nvmu/I");
  trout->Branch("nve",&vis->nve,"nve/I");
  trout->Branch("nvgam",&vis->nvgam,"nvgam/I");
  trout->Branch("nvpip",&vis->nvpip,"nvpip/I");
  trout->Branch("nvpi0",&vis->nvpi0,"nvpi0/I");
  trout->Branch("nvp",&vis->nvp,"nvp/I");
  trout->Branch("nvk",&vis->nvk,"nvk/I");
  trout->Branch("visbrightness",vis->visbrightness,"visbrightness[100]/D");
  trout->Branch("viswall",vis->viswall,"viswall[100]/D");
  trout->Branch("vistowall",vis->vistowall,"vistowall[100]/D");
  trout->Branch("vismrwall1",&vis->vismrwall1,"vismrwall1/D");
  trout->Branch("vismrwall2",&vis->vismrwall2,"vismrwall2/D");
  trout->Branch("vismrtowall1",&vis->vismrtowall1,"vismrtowall1/D");
  trout->Branch("vismrtowall2",&vis->vismrtowall2,"vismrtowall2/D");
  trout->Branch("vismrwallmin",&vis->vismrwallmin,"vismrwallmin/D");
  trout->Branch("vismrtowallmin",&vis->vismrtowallmin,"vismrtowallmin/D");
  trout->Branch("nvisscnd",&vis->nvisscnd,"nvisscnd/I");
  trout->Branch("vismrbrightness",&vis->vismrbrightness,"vismrbrightness/D");
  trout->Branch("vismrpid1",&vis->vismrpid1,"vismrpid1/I");
  trout->Branch("vismrpid2",&vis->vismrpid2,"vismrpid2/I");
  trout->Branch("vismrt1",&vis->vismrt1,"vismrt1/D");
  trout->Branch("vismrt2",&vis->vismrt2,"vismrt2/D");
  trout->Branch("vismrtype1",&vis->vismrtype1,"vismrtype1/I");
  trout->Branch("vismrtype2",&vis->vismrtype2,"vismrtype2/I");
  trout->Branch("vistime",vis->vistime,"vistime[100]/D");
  trout->Branch("vispid",vis->vispid,"vispid[100]/I");
  trout->Branch("visscndpid",vis->visscndpid,"visscndpid[100]/I");
  trout->Branch("visscndparentid",vis->visscndparentid,"visscndparentid[100]/I");
  // fitqun
  trout->Branch("fqwall",&wall,"fqwall/F");
  trout->Branch("fqtowall",&towall,"fqtowall/F");
  trout->Branch("fq1rwall",fq1rwall,"fq1rwall[10][7]/F");
  trout->Branch("fq1rtowall",fq1rtowall,"fq1rtowall[10][7]/F");
  trout->Branch("towallv",towallv,"towallv[50]");
  trout->Branch("wallv2",&wallv2,"wallv2");
  trout->Branch("evtweight",&evtweight,"evtweight/F");
  trout->Branch("best2RID",&best2RID,"best2RID/I");
  trout->Branch("fq1rperim",fq1rperim,"fq1rperim[10][7]/F");
  trout->Branch("fq1rmincone",fq1rmincone,"fq1rmincone[10][7]/F");

  return;
}


/////////////////////////////
//empty constructor
preProcess::preProcess(){
  nFiles=0;
  useWeights=0;
}


//////////////////////////////////////////
//read in parameters and run preprocessing!
void preProcess::runPreProcessing(){

   //process the files
  TString datatag = nameTag.Data();
  datatag.Append("_ppdata");
  TString mctag = nameTag.Data();
  mctag.Append("_ppmc");
  nameTag = mctag.Data();
  processAllFiles(chmc);
  nameTag = datatag.Data();
  processAllFiles(chdat); 

  cout<<"preProcess: Complete!"<<endl;

  //////////////////////////
  return;
}

#endif
