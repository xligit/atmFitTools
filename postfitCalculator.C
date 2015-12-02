//#include "TTree.h"
//#include "TFile.h"
//#include "TString.h"
//#include "TH1F.h"

//#include "shared.h"

#include "atmFitPars.C"
#include "histoManager.C"
//#include "fQreader.C"
#include "splineFactory.C"
#include "mcmcReader.C"

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
  float pars[1000]; 
  float evtweight;
  
  TString parFileName;
  //histogrms
  TH1F* hEnu;
 

  ////////////////////////////////////////////////////////
  //methods
  void  setMCTree(TTree* tr);
  float getEvtWeight(); //< returns the weight given the current parameters
  void  setParsFromMCMC(int istep); //< sets parameters from mcmmc cloud
  void fillHistos(); //fills all histograms
  int makeSelection(int iselect); //select events based on various attributes;
};

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

void postfitCalculator::fillHistos(){
  hEnu->Fill(mcreader->attribute[3],getEvtWeight());
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
float  postfitCalculator::getEvtWeight(){

  /////////////////////////////////////////////
  //get event weights
  float ww = 1.0;
  for (int isyspar=0;isyspar<fitpars->nSysPars;isyspar++){
    ww*=sfact->getEvtWeight(mcreader,isyspar,fitpars->sysPar[isyspar]);
  } 
  evtweight = ww;

  ////////////////////////////////////////////////////
  //modify attributes
  for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
    mcreader->attribute[iatt]*=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0]; //multiply by smear parameter
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
  float xmin =  0.;
  float xmax = 1000.;
  hEnu = new TH1F("enu","enu",nbins,xmin,xmax);

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

