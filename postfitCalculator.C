//#include "TTree.h"
//#include "TFile.h"
//#include "TString.h"
//#include "TH1D.h"

//#include "shared.h"

#include "atmFitPars.C"
#include "histoManager.C"
//#include "fQreader.C"
#include "splineFactory.C"
#include "mcmcReader.C"
#include "THStack.h"
#include "sharedPars.C"


#define NCLASSMAX 10

//This class aims to take the MCMC output and generate uncertainties in number of events 
class postfitCalculator{
  public:
  ///////////////////////////////////////////////////////////////////////////////
  //constructor
  //initialize using the file name that the contains the mcmc tree to be analyzed
  //the optional parameter is the name of the shared paremeter file for creating the
  //atmFitPars object.
  postfitCalculator(const char* mcmcfilename,const char* parfile);
  postfitCalculator(const char* parfile); //< construct from parameter file
  void init(); //< called by constructor to initialize histograms, trees, etc

  /////////////////////////////////////////////////////////////
  //internal variables
  TTree* mcmcpath; //< points to tree containing the mcmc steps
  TTree* mctree;  //< points to the tree containg the mc events
  TTree* datatree; //< for data events
  fQreader* datareader; //< points t
  mcmcReader* path; //< points to a step in mcmcpath
  fQreader* mcreader; //< points to an event in mctree
  atmFitPars* fitpars; //< points to the current fit parameters object
  splineFactory* sfact; //< empty spline factory to call getEvtWeight. 
  TFile* outFile; //< outputfile
  int NMCEvents; //< limits # of MC events to use to speed up histo filling
  double pars[1000]; 
  double evtweight;
  int evtclass; //< code for event class 
  TString parFileName;  //< name of parameter file
  //for "nsk" tree
  double nevents[10]; //< number of events for each event class
  double neventsdef[10]; //< number of events for each class in default mc
  double neventstot[10]; //< total number of events for each event class
  TTree* nsktree; //< tree to hold # of events for each MCMC point
  int selectiontype; //< code for which selection to use
  int NMCMCPts; //< number of points to sample from
  int MCMCBurnIn; //< number of burn-in points
  sharedPars* runPars;
  int currentMCMCPoint;
  int currentMCEvent;

  //histogrms of interest
  TH1D* hEnuMu[10][200]; //<energy assuming muon
//  TH1D* hPIDemuSE2[10][200]; //< e/mu PID of second subev
  TH1D* hPIDemu[10][200]; //<e/mu PID
  TH1D* hNSK[10]; //< numbers of events after cuts
  

  ////////////////////////////////////////////////////////
  //methods
  // modfies the values of the current MC event using parameters from current MCMC point
  void  modifyCurrentEvent();
  void  setMCTree(TTree* tr);
  void  setSelectionType(int itype){selectiontype=itype;}
  double getEvtWeight(); //< returns the weight given the current parameters
  void  modifyAttributes();
  void  unmodifyAttributes();
  void  setParsFromMCMC(int istep); //< sets parameters from mcmmc cloud
  void fillHistos(int iclass,int islot); //< fills all histograms 
  int makeSelection(int iselect); //< select events based on various attributes;
  int getEvtClass();
  void drawBreakdown(int ihisto,int islot);
  void makeHistos(int iselection,int islot);
  void findEvtClasses(); //< fills array of total numbers of each event class
  void fillNskTree(int npts,int nburn=1000); //< step through MCMC cloud and fill NSK tree at each point
  void makeNSKTable();
  void showAllHistos(int ihisto,int iclass);
  void MCMCLooperTemplate();
  void cosmicPostFitAnalysis();
//  void runPostFit();
};



void  postfitCalculator::cosmicPostFitAnalysis(){


  //get number of total steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-MCMCBurnIn;
  if (nsteps<(NMCMCPts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  //loop over points (uniform sampling)
  int dstep = nsteps/NMCMCPts;
  int istep = MCMCBurnIn;
  int samppts[10000];

  //get array of points to sample from
  for (int i=0;i<NMCMCPts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

    //find the number of MC events to use
  if ((NMCEvents<=0)||(NMCEvents>mctree->GetEntries())) NMCEvents = mctree->GetEntries(); 

  int NmuID[100];
  int NmuMisID[100];
  int NeID[100];
  int NeMisID[100];
 
  //set number of mis ied events to zer
  for (int i=0;i<100;i++){
    NmuID[i]=0.;
    NmuMisID[i]=0.;
    NeID[i]=0.;
    NeMisID[i]=0.;
  }
 
  //loop over number of MCMC points
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    currentMCMCPoint = samppts[ipt];
    //loop over MC events
    for (int ievt=0;ievt<NMCEvents;ievt++){
      currentMCEvent = ievt;
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      if (ipt!=0) modifyCurrentEvent();//< all attributes are modified and eventWeight is calculated
      //feel free to do things here, like fill histograms or something  
      hPIDemu[0][ipt]->Fill(mcreader->attribute[0]);
      hPIDemu[1][ipt]->Fill(mcreader->attribute[1]);
      if (makeSelection(1)){
        NmuID[ipt]++;
      }
      else{
        NmuMisID[ipt]++;
      }
      if (makeSelection(2)){
        NeID[ipt]++;
      }
      else{
        NeMisID[ipt]++;
      }
    } 
  }

  double muratemean=0.;
  double eratemean=0.;
  double muratevar=0.;
  double eratevar=0.;
  for (int jpt=0;jpt<NMCMCPts;jpt++){
    float muRate = (float)NmuMisID[jpt]/((float)NmuMisID[jpt]+(float)NmuID[jpt]);
    float eRate = (float)NeMisID[jpt]/((float)NeMisID[jpt]+(float)NeID[jpt]);
    cout<<" NmuID  "<<NmuID[jpt]<<" NmuMISID  "<<NmuMisID[jpt]<<" NeID  "<<NeID[jpt]<<" NeMisID  "
    <<NeMisID[jpt]<<" muRate  "<<muRate<<" eRate  "<<eRate<<endl;
    if (jpt!=0) muratemean+=muRate;
    if (jpt!=0) eratemean+=eRate;
  }

  muratemean/=(float)(NMCMCPts-1);
  eratemean/=(float)(NMCMCPts-1);

  for (int jpt=0;jpt<NMCMCPts;jpt++){
    float muRate = (float)NmuMisID[jpt]/((float)NmuMisID[jpt]+(float)NmuID[jpt]);
    float eRate = (float)NeMisID[jpt]/((float)NeMisID[jpt]+(float)NeID[jpt]);
    if (jpt!=0) muratevar+=(muRate-muratemean)*(muRate-muratemean);
    if (jpt!=0) eratevar+=(eRate-eratemean)*(eRate-eratemean);
  }

 
  muratevar=TMath::Sqrt(muratevar);
  eratevar=TMath::Sqrt(muratevar);
  muratevar/=(float)(NMCMCPts-2);
  eratevar/=(float)(NMCMCPts-2);



  cout<<"avg muon misID rate: "<<muratemean<<endl;
  cout<<"muon misID var: "<<muratevar<<endl;
  
  cout<<"avg electron misID rate: "<<eratemean<<endl;
  cout<<"electron misID var: "<<eratevar<<endl;


 
  return;

}

//void postfitCalculator::runPostFit(){

 // return;
//}

postfitCalculator::postfitCalculator(const char* parfile){
  
  //read in parameters
  runPars = new sharedPars(parfile);
  runPars->readParsFromFile();

  NMCMCPts = runPars->NMCMCPts;  
  MCMCBurnIn = runPars->MCMCBurnIn;
  NMCEvents = runPars->NMCEvents;
  //set parameter file name
  parFileName = parfile;

  //get mcmc path tree
  TString mcmcfilename = runPars->MCMCFile;
  TFile *mcmcfile = new TFile(mcmcfilename.Data());
  mcmcpath = (TTree*)mcmcfile->Get("MCMCpath");
  path = new mcmcReader(mcmcpath);


  //setup atmFitpars
  fitpars = new atmFitPars(parfile); 

  //initialize histograms
  init();    


}

/////////////////////////////////////////
//draws all modified histograms on same pad
void postfitCalculator::showAllHistos(int ihisto, int iclass){
  
  //loop over all histograms
  hEnuMu[iclass][0]->SetLineColor(kRed);
  if (ihisto==0)  hEnuMu[iclass][0]->Draw();
  if (ihisto==1)  hPIDemu[iclass][0]->Draw();
  for (int i=1;i<NMCMCPts;i++){
    if (ihisto==0){
      hEnuMu[iclass][i]->Draw("same");
    }
    if (ihisto==1){
      hPIDemu[iclass][i]->Draw("same");
    }
  }
  
  return;
}

/////////////////////////////////////////
//Fill table showing event numbers
void postfitCalculator::makeNSKTable(){

  //arrays to be filled
  double NSK[NCLASSMAX];
  double NSKRMS[NCLASSMAX];

  nsktree->SetBranchAddress("nevents",nevents);
  nsktree->SetBranchAddress("neventstot",neventstot);
  nsktree->SetBranchAddress("neventsdef",neventsdef);
  
  //set initial to zero
  for (int j=0;j<NCLASSMAX;j++){
    NSK[j]=0.;
    NSKRMS[j]=0.;
  }

  //loop over points to get mean
//  for (int i=0;i<nsktree->GetEntries();i++){
  for (int i=0;i<10;i++){
    nsktree->GetEntry(i);
    for (int j=0;j<NCLASSMAX;j++){
      NSK[j]+=nevents[j];
    }
  }
  for (int j=0;j<NCLASSMAX;j++){
    if (NSK[j]>=1.0){
      cout<<"NSK: "<<NSK[j]<<endl;
      NSK[j]/=10.; //< normalize for average # of events
    }
  }

  //loop over points to get RMS
//  for (int i=0;i<nsktree->GetEntries();i++){
  for (int i=0;i<10;i++){
    nsktree->GetEntry(i);
    for (int j=0;j<NCLASSMAX;j++){
      if (nevents[j]>0.1){
        NSKRMS[j] += ((nevents[j]-NSK[j])*(nevents[j]-NSK[j])); //< add square deviation
      }
    }
  }
  for (int j=0;j<NCLASSMAX;j++){
    if (NSKRMS[j]>=1.0){
      cout<<"NSKRMS: "<<NSKRMS[j]<<endl;      
      NSKRMS[j]/=10.; //< normalize 
      NSKRMS[j] = TMath::Sqrt( NSKRMS[j] ); //< get root mean square
    }
  }
  
  cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  for (int j=0;j<NCLASSMAX;j++){
    cout<<"Class: "<<j
    <<" Ndef: "<<neventsdef[j]
    <<" N: "<<NSK[j]
    <<" NRMS: "<<NSKRMS[j]
    <<" n: "<<(NSK[j]-neventsdef[j])/neventsdef[j]
    <<" nRMS: "<<NSKRMS[j]/neventsdef[j]
    <<endl;
  }
  cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;

  return;
}

///////////////////////////////////////////////////
//undo attribute modification
void postfitCalculator::unmodifyAttributes(){

  ////////////////////////////////////////////////////
  //un-modify attributes
  for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
    //mcreader->attribute[iatt]/=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0]; //divide by smear parameter
    mcreader->attribute[iatt]-=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1]; //subtract bias parameter
  }

  //////////
  return;
}


///////////////////////////////////////////////////
//changes attributes of event to the modified quantities
void postfitCalculator::modifyAttributes(){
//  cout<<"N "<<fitpars->nAttributes<<endl;
  ////////////////////////////////////////////////////
  //modify attributes
  for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
//    cout<<"attribute mod: "<<fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1]<<endl;
    //mcreader->attribute[iatt]*=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0]; //multiply by smear parameter
    mcreader->attribute[iatt]+=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1]; //add bias parameter
  }

  //////////
  return;
}

///////////////////////////////////////////////
// modfies the values of the current MC event using parameters from current MCMC point
void postfitCalculator::modifyCurrentEvent(){
//  cout<<"test"<<endl;
  setParsFromMCMC(currentMCMCPoint); //< sets parameters from current step of MCMC path

  mctree->GetEntry(currentMCEvent); //< read in default MC event

  getEvtWeight(); //< gets the event weight form the current flux and xsec parameters

  modifyAttributes(); //< modifies each attribute according to the current MCMC point 
  
  return;
}

/////////////////////////////////////////////////////////////
//template for looping over MCMC points and modifying the MC
void postfitCalculator::MCMCLooperTemplate(){

  //get number of total steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-MCMCBurnIn;
  if (nsteps<(NMCMCPts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  //loop over points (uniform sampling)
  int dstep = nsteps/NMCMCPts;
  int istep = MCMCBurnIn;
  int samppts[10000];

  //get array of points to sample from
  for (int i=0;i<NMCMCPts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

  //find the number of MC events to use
  if ((NMCEvents<=0)||(NMCEvents>mctree->GetEntries())) NMCEvents = mctree->GetEntries(); 

  //loop over number of MCMC points
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    currentMCMCPoint = samppts[ipt];
    //loop over MC events
    for (int ievt=0;ievt<NMCEvents;ievt++){
      currentMCEvent = ievt;
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      modifyCurrentEvent();//< all attributes are modified and eventWeight is calculated

      //feel free to do things here, like fill histograms or something  
 
    } 
  }
 
  return;

}



/////////////////////////////////////////////////
//Fill the NSK tree at various steps in MCMC path
void postfitCalculator::fillNskTree(int npts,int nburn){

  NMCMCPts = npts;
  
  //get number of steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-nburn;
  if (nsteps<(npts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  //loop over points (uniform sampling)
  int dstep = nsteps/npts;
  int istep = nburn;
  int samppts[10000];
  //get array of sampled points
  for (int i=0;i<npts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

  //get total numbers of events
  findEvtClasses();


  //loop over number of points
  for (int ipt=0;ipt<npts;ipt++){
    setParsFromMCMC(samppts[ipt]);//< fill atmFitPars from parameters at this step
    for (int ievt=0;ievt<NMCEvents;ievt++){
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      mctree->GetEntry(ievt);
      evtclass = getEvtClass(); //< get event class
      evtweight = getEvtWeight(); //< gets weight for event and modifies attributes   
      modifyAttributes(); //<  modify the attributes according the the mcmc parameters
      if (makeSelection(selectiontype)){
        fillHistos(evtclass,ipt); //< fills the histograms
        nevents[evtclass]++; //< increment the number of events of this class passing the cut
      }
    } 
 
    //fill nsk
    nsktree->Fill();

    //reset nevents
    for (int iclass=0;iclass<NCLASSMAX;iclass++){
      nevents[iclass]=0.;
    }
  }

  //write out output 
  nsktree->Write();
  //outFile->Close();

  ////////////////
  return; 
}

//////////////////////////////////////////////
//Find total numbers of events of each class
void postfitCalculator::findEvtClasses(){

  //first make sure each event class has zero events
  for (int iclass=0;iclass<NCLASSMAX;iclass++){
    neventstot[iclass]=0.;
    nevents[iclass]=0.;  
  }

  //now loop through all MC events to see hou many events of each
  //class we have
  cout<<"Finding total # of events of each class...."<<endl;
  for (int ievt=0;ievt<NMCEvents;ievt++){
    mctree->GetEntry(ievt);
    evtclass=getEvtClass(); //< fills "evtclass" variable
    neventstot[evtclass]++; //<increment the total number of events of this class
    if (makeSelection(selectiontype)){
      neventsdef[evtclass]++; //< increment the number of events of this class passing the cut
    }
  }
 
 
  /////////
  return;
}

//////////////////////////////////////////////
//Draw all event classes for a histogram
void postfitCalculator::drawBreakdown(int ihisto,int islot){

  ///////////////////
  //pid histograms
  if (ihisto==0){
    hPIDemu[0][islot]->SetFillColor(kBlue);
    hPIDemu[1][islot]->SetFillColor(kRed);
    hPIDemu[2][islot]->SetFillColor(kCyan);
    hPIDemu[3][islot]->SetFillColor(kOrange);
    hPIDemu[4][islot]->SetFillColor(kBlack); 
  }

  ///////////////////
  //Enu (mu) histograms
  if (ihisto==1){
    hEnuMu[0][islot]->SetFillColor(kBlue);
    hEnuMu[1][islot]->SetFillColor(kRed);
    hEnuMu[2][islot]->SetFillColor(kCyan);
    hEnuMu[3][islot]->SetFillColor(kOrange);
    hEnuMu[4][islot]->SetFillColor(kBlack);    
  }

  ////////////
  return;
}

/////////////////////////////////////////////////////////////
//returns an integer corresponding to the true event class
int postfitCalculator::getEvtClass(){
  if ((TMath::Abs(mcreader->mode)==1)&&(TMath::Abs(mcreader->ipnu[0])==12)) return 0; //< CCQE nu-e
  if ((TMath::Abs(mcreader->mode)==1)&&(TMath::Abs(mcreader->ipnu[0])==14)) return 1; //< CCQE nu-mu
  if ((TMath::Abs(mcreader->mode)>1)&&
      (TMath::Abs(mcreader->mode)<30)&&
      (TMath::Abs(mcreader->ipnu[0])==12)) return 2; //< CCnQE nu-e
  if ((TMath::Abs(mcreader->mode)>1)&&
      (TMath::Abs(mcreader->mode)<30)&&
      (TMath::Abs(mcreader->ipnu[0])==14)) return 3; //< CCnQE nu-mu
  if ((TMath::Abs(mcreader->mode)>30)) return 4; //< NC

  //////////////
  return -1;
}

void postfitCalculator::makeHistos(int iselect, int islot){
  
  //loop over mc events
  for (int ievt=0;ievt<mctree->GetEntries();ievt++){
    if ((ievt%1000)==0) cout<<ievt<<endl;
    mctree->GetEvent(ievt); //< load tree info
    evtweight = getEvtWeight(); //< gets weight for event and modifies attributes   
    modifyAttributes();
    evtclass = getEvtClass(); //< get event class
    if (makeSelection(iselect)){
      fillHistos(evtclass,islot); //< fills the histograms
      nevents[evtclass]++; //< increment the number of events of this class passing the cut
    }
  }

  //fill nsk tree for this point  
//  nsktree->Fill();

  ////////////
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//make various event selections, indexed by "iselect"
//returns 1 if all cuts are passed, 0 otherwise
int postfitCalculator::makeSelection(int iselect){

  ///////////////////////////
  //mock 1R numu calculation
  if (iselect==0){
    //Evis cut
    if (mcreader->attribute[2]<100.) return 0;
    //FC cut
    if (mcreader->nhitac>16) return 0;
    //PID cut
    if (mcreader->attribute[0]>0.) return 0;
    //nring cut
    if (mcreader->fqmrnring[0]!=1) return 0;
  }


  /////////////////////////////
  //cosmic mu mock selection
  if (iselect==1){
    //PID cut
    if ((-1.*mcreader->attribute[0])<0.2*mcreader->fq1rmom[0][1]) return 0;
    //nring cut
  //  if (mcreader->attribute[2]>0) return 0;
    //evis cut
    if (mcreader->fq1rmom[0][1]<30) return 0;
  }

  /////////////////////////////
  //cosmic decay e mock selection
  if (iselect==2){
    //PID cut
    if ((-1.*mcreader->attribute[1])>0.2*mcreader->fq1rmom[1][1]) return 0;
    //nring cut
  //  if (mcreader->attribute[2]>0) return 0;
    //evis cut
    if (mcreader->fq1rmom[0][1]<30) return 0;
  }

  //event has passed cuts
  return 1;
}


void postfitCalculator::fillHistos(int iclass,int islot){
  hEnuMu[iclass][islot]->Fill(mcreader->attribute[3],evtweight);
  hPIDemu[iclass][islot]->Fill(mcreader->attribute[0],evtweight);

}


///////////////////////////////////////////////////////////////////////////////////
//fills the atmFitPars object with the parameters from step "istep" of the MCMC
void  postfitCalculator::setParsFromMCMC(int istep){
  mcmcpath->GetEntry(istep); //< fills mcmcpath reader
  for (int ipar=0;ipar<path->npars;ipar++){
    fitpars->setParameter(ipar,path->par[ipar]);
  } 
  return;
}


void postfitCalculator::setMCTree(TTree* tr){
  mctree = tr;
  mcreader = new fQreader(mctree);
  mctree->GetEntry(0);
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//get event weight according to the current flux and xsec parameters
double  postfitCalculator::getEvtWeight(){

  /////////////////////////////////////////////
  //get event weights
  double ww = 1.0;
  //loop over all flux and xsec parameters
  for (int isyspar=0;isyspar<fitpars->nSysPars;isyspar++){
    ww*=sfact->getEvtWeight(mcreader,isyspar,fitpars->sysPar[isyspar]);
  } 
  evtweight = ww;

  ///////////////
  return ww;
}


///////////////////////////////////////////////////////////////////////////////
//inialize some values
void postfitCalculator::init(){

  //////////////////////
  //make tree for numbers of events
  outFile = new TFile("nsk.root","recreate");
  nsktree = new TTree("nsk","nsk");
  nsktree->Branch("par",path->par,"par[100]/D");
  nsktree->Branch("nevents",nevents,"nevents[10]/D");
  nsktree->Branch("neventsdef",neventsdef,"neventsdef[10]/D");
  nsktree->Branch("neventstot",neventstot,"neventstot[10]/D");

  //////////////////////  
  //setup histograms
  int nbins=25;
  double xmin =  0.;
  double xmax = 5000.;
  double xminpid = -5000;
  double xmaxpid = 5000;
  for (int iclass=0;iclass<10;iclass++)
  {
    for (int i=0;i<200;i++){
      hEnuMu[iclass][i] = new TH1D(Form("enu_class%d_%d",iclass,i),Form("enu_class%d_%d",iclass,i) ,nbins,xmin,xmax);
      hPIDemu[iclass][i] = new TH1D(Form("PIDemu_class%d_%d",iclass,i),Form("PIDemu_class%d_%d",iclass,i) ,nbins,xminpid,xmaxpid);
    }
  }

  ////////////////////////////
  //setup atm fit pars
//  fitpars = new atmFitPars(parFileName.Data());
//  fitpars->initPars("tn186");

  ////////////////////////////
  //setup spline factory
  sfact = new splineFactory();  

  ////////////////////////////
  //other defaults
  NMCEvents=0;
  selectiontype=0;   
  return;
}


postfitCalculator::postfitCalculator(const char* mcmcfilename,const char* parfile){
  
  //set parameter file name
  parFileName = parfile;
  //get mcmc path
  TFile *mcmcfile = new TFile(mcmcfilename);
  mcmcpath = (TTree*)mcmcfile->Get("MCMCpath");
  path = new mcmcReader(mcmcpath);
  init();    
}

