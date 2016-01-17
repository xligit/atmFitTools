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

//This class aims to take the MCMC output and generate uncertainties in number of events 
class postfitCalculator{
  public:
  ///////////////////////////////////////////////////////////////////////////////
  //constructor
  //initialize using the file name that the contains the mcmc tree to be analyzed
  //the optional parameter is the name of the shared paremeter file for creating the
  //atmFitPars object.
  postfitCalculator(const char* mcmcfilename,const char* parfile="");
  void init(); //< called by constructor to initialize histograms, trees, etc

  /////////////////////////////////////////////////////////////
  //internal variables
  TTree* mcmcpath; //< points to tree containing the mcmc steps
  TTree* mctree;  //< points to the tree containg the mc events
  mcmcReader* path; //< points to data in mcmcpath
  fQreader* mcreader; //< points to data in mctree
  atmFitPars* fitpars; //< points to the current fit parameters object
  splineFactory* sfact; //< empty spline factory to call getEvtWeight. 
  double pars[1000]; 
  double evtweight;
  int evtclass; //< code for event class 
  TString parFileName;

  //histogrms
  TH1D* hEnuMu[10][200]; //<energy assuming muon
  TH1D* hPIDemu[10][200]; //<e/mu PID

  ////////////////////////////////////////////////////////
  //methods
  void  setMCTree(TTree* tr);
  double getEvtWeight(); //< returns the weight given the current parameters
  void  setParsFromMCMC(int istep); //< sets parameters from mcmmc cloud
  void fillHistos(int iclass,int islot); //fills all histograms with data pointed to by mcreader
  int makeSelection(int iselect); //select events based on various attributes;
  int getEvtClass();
  void drawBreakdown(int ihisto,int islot);
  void makeHistos(int iselection,int islot);
};

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
  return -1;
}

void postfitCalculator::makeHistos(int iselect, int islot){
  
  //loop over mc events
  for (int ievt=0;ievt<mctree->GetEntries();ievt++){
    if ((ievt%1000)==0) cout<<ievt<<endl;
    mctree->GetEvent(ievt); //< load tree info
    evtweight = getEvtWeight(); //< gets weight for event and modifies attributes   
    evtclass = getEvtClass(); //< get event class
    if (makeSelection(iselect)) fillHistos(evtclass,islot); //<fills the histograms
  }

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
//get event weights and modify attributes according to parameters
double  postfitCalculator::getEvtWeight(){

  /////////////////////////////////////////////
  //get event weights
  double ww = 1.0;
  for (int isyspar=0;isyspar<fitpars->nSysPars;isyspar++){
    ww*=sfact->getEvtWeight(mcreader,isyspar,fitpars->sysPar[isyspar]);
  } 
  evtweight = ww;

  ////////////////////////////////////////////////////
  //modify attributes
  for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
//    mcreader->attribute[iatt]*=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0]; //multiply by smear parameter
    mcreader->attribute[iatt]+=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1]; //add bias parameter
  }
 
  return ww;
}


///////////////////////////////////////////////////////////////////////////////
//inialize some values
void postfitCalculator::init(){

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
  fitpars = new atmFitPars(parFileName.Data());
  fitpars->initPars("tn186");

  ////////////////////////////
  //setup spline factory
  sfact = new splineFactory();  
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

