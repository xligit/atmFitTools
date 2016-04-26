#ifndef ATMFITPARS_C
#define ATMFITPARS_C

#include "shared.h"
#include "TRandom2.h"
#include "sharedPars.C"

using namespace std;

//class containing all atmospheric fit pars
class atmFitPars{
  public:
  
  ///////////////////////////////////////////////////////////////////
  //constructors
  atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst=0);
  atmFitPars(int isamp, int ibin, int icomp, int iatt, const char* systype); 
  // use this constructor
  atmFitPars(const char* parfile); //constructs from parameter file

  ///////////////////////////////////////////////////////////////
  //numbers of various parametrs
  int nSamples;
  int nBins;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nTotPars;
  int flgUseNormPars;
  sharedPars* runpars;

  ///////////////////////////////////////////////////////////////////
  //parameter values
  double histoNorm[NSAMPMAX][NBINMAX];
  double histoPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncLo[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  double sysPar[NSYSPARMAX];
  double sysParDefault[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  double pars[4000];
  double parUnc[4000];
  int   fixPar[4000]; //< array of fix flags for parameters
  double bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2]; //< stores 1D array position for bias/smear pars
  int   sysParIndex[NSYSPARMAX]; //< stores 1D array position for systematic pars
  int   normParIndex[NSAMPMAX][NBINMAX]; //< stores 1D array position for normalization pars

  double norm;  

  //////////////////////////////////////////////////////////////
  //methods
  void setNorm(double x){norm=x;}
  void initPars(const char* systype=""); //< sets parameters to initial values
  int getParIndex(int ibin, int icomp, int iatt, int imod){return parIndex[ibin][icomp][iatt][imod];}
  double getParameter(int ipar){return pars[ipar];}
  double getHistoParameter(int ibin, int icomp, int iatt, int imod);
  double getSysParameter(int isys);
  double getNormParameter(int isamp, int ibin);
  void setParameter(int ipar, double value);
  void setSysParameter(int ipar, double value);
  void setParameter(int ibin, int icomp, int iatt, int imod, double value); 
  void setSysParUnc(int isys,double value){sysParUnc[isys]=value;}
  void fixParameter(int ipar);
  void fixParameter(int ibin,int icomp,int iatt, int imod);
  void fixAllSmearPars(int isfixed=1);
  void setRandSysPar(); //sets systematic parameters to random values
  int  checkFixFlg(int ibin,int icomp,int iatt, int imod);
  void resetDefaults();
  void printParValues();
  int binOfPar[4000];
  int compOfPar[4000];
  int attOfPar[4000];
  int typeOfPar[4000];

  //saving and reading pars
  void savePars(const char* filename);
  void readPars(const char* filename);
  void printPars();
};

//////////////////////////////////////////////////
// get specific parameters
double atmFitPars::getHistoParameter(int ibin, int icomp, int iatt, int imod){
  int theindex = parIndex[ibin][icomp][iatt][imod];
  return pars[theindex];
}

//////////////////////////////////////////////////
// get specific parameters
double atmFitPars::getNormParameter(int isamp, int ibin){
  int theindex = normParIndex[isamp][ibin];
  return pars[theindex];
}

//////////////////////////////////////////////////
// get specific parameters
double atmFitPars::getSysParameter(int isys){
  int theindex = sysParIndex[isys];
  return pars[theindex];
}

//////////////////////////////////////////////////
//set parameters back to defaults
void atmFitPars::resetDefaults(){

//  initPars();

  
  //initialize histogram pars
  int index = 0; //< running 1D index
  for (int ibin=0;ibin<nBins;ibin++){
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        histoPar[ibin][icomp][iatt][0]=1.0;
        pars[index]=1.0;
        parUnc[index]=0.005; //rough estimate of uncertainty
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=0;
        fixPar[index]=0;
        parIndex[ibin][icomp][iatt][0]=index;  
        index++;
        histoPar[ibin][icomp][iatt][1]=0.0;
        parIndex[ibin][icomp][iatt][1]=index;
        pars[index]=0.0;
        parUnc[index]=1.0;  //rough estimate of uncertainty
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=1;
        fixPar[index]=0;
        index++;
      }
    }
  }

  //initialize systematic error parameters
  for (int isyst=0;isyst<nSysPars;isyst++){
    pars[index]=sysParDefault[isyst];
    index++;
    sysPar[isyst]=sysParDefault[isyst];
  }
  

  return;
}


//////////////////////////////////////////////////
//Set only bias and systematic pars to float in fit
void atmFitPars::fixAllSmearPars(int isfixed){
  //fix all smear parameters
  for (int ibin=0;ibin<nBins;ibin++){
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        int index = parIndex[ibin][icomp][iatt][0];
        fixPar[index] = isfixed;
      }
    }
  }
  return;
}

/////////////////////////////////////////////////////
//construct from parameter file
atmFitPars::atmFitPars(const char* parfilename){

  //fill shared parameters from file
  cout<<"atmFitPars: reading parameter file: "<<parfilename<<endl;
  runpars = new sharedPars(parfilename);
  runpars->readParsFromFile();
  nSamples = runpars->nSamples;
  cout<<"atmFitPars:  nSamples: "<<nSamples<<endl;
  nComponents = runpars->nComponents;
  cout<<"atmFitPars:  nComponents: "<<nComponents<<endl;
  nBins = runpars->nFVBins;
  cout<<"atmFitPars:  nBins: "<<nBins<<endl;
  nSysPars = runpars->nSysPars;
  cout<<"atmFitPars:  nSysPars: "<<nSysPars<<endl;
  nAttributes = runpars->nAttributes;
  cout<<"atmFitPars:  nAttributes: "<<nAttributes<<endl;
  flgUseNormPars = runpars->flgUseNormPars;
  if (flgUseNormPars) cout<<"atmFitPars: Using normalization pars"<<endl;
  TString sysType = runpars->sysParType;

  //fill all initial parameter values and count number of pars
  initPars(sysType.Data());

}

void atmFitPars::printPars(){
  cout<<"$$$ CURRENT PARAMETER VALUES $$$"<<endl;
  for (int i=0;i<nTotPars;i++){
    cout<<"PAR: "<<i<<" = "<<pars[i]<<" +/- "<<parUnc[i]<<endl;
  }
  cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  return;
}


//read in all parameters from a file made from a call to "savePars"
void atmFitPars::readPars(const char* filename){
  //open parameter file
  TFile* fpars = new TFile(filename);
  //get parameter tree
  TTree* parTree = (TTree*)fpars->Get("parTree");
  parTree->SetBranchAddress("nTotPars",&nTotPars);
  parTree->SetBranchAddress("nSysPars",&nSysPars);
  parTree->SetBranchAddress("pars",pars);
  parTree->SetBranchAddress("parUnc",parUnc);
  parTree->GetEntry(0);
  //set parameter explicitly to make sure arrays are filled as well
  for (int ipar=0;ipar<nTotPars;ipar++){
    cout<<"setting parameter # "<<ipar<<" to "<<pars[ipar]<<endl;
    setParameter(ipar,pars[ipar]);
  }
  return;
}

//save current parameters to a file 
void atmFitPars::savePars(const char* filename){
  //create output file
  TFile* fout = new TFile(filename,"RECREATE");
  //create tree to hold parameter values
  TTree* parTree = new TTree("parTree","parTree");
  parTree->Branch("pars",pars,"pars[4000]/D");
  parTree->Branch("parUnc",parUnc,"parUnc[4000]/D");
  parTree->Branch("nTotPars",&nTotPars,"nTotPars/I");
  parTree->Branch("nSysPars",&nSysPars,"nSysPars/I");
  parTree->Fill();
  fout->Write();
  return;
}

TRandom2* randy2 = new TRandom2();

void atmFitPars::setRandSysPar(){
 // TRandom2* randy = new TRandom2();
  double parval;
  for (int i=0;i<nSysPars;i++){
     parval = randy2->Gaus(sysPar[i],(sysParUnc[i]/2.));
     if (parval<0) parval=0.;
     cout<<"par "<<i<<" is "<<parval<<endl;
     setSysParameter(i,parval);
  }
  return;
}

void atmFitPars::printParValues(){
  for (int ipar=0;ipar<nTotPars;ipar++){
    cout<<"par "<<ipar<<": "<<pars[ipar]<<endl;
  }
}

int atmFitPars::checkFixFlg(int ibin, int icomp, int iatt, int imod){
  return fixPar[parIndex[ibin][icomp][iatt][imod]];
}

void atmFitPars::fixParameter(int ibin,int icomp,int iatt, int imod){
  fixPar[getParIndex(ibin,icomp,iatt,imod)] = 1;
  return;
}

void atmFitPars::fixParameter(int ipar){
  fixPar[ipar] = 1;
  return;
}

void atmFitPars::setSysParameter(int ipar, double value){
  sysPar[ipar]=value;
  pars[nTotPars-nSysPars+ipar]=value;
  return;
}

void atmFitPars::setParameter(int ibin, int icomp, int iatt, int itype, double value){
  histoPar[ibin][icomp][iatt][itype] = value;
  pars[parIndex[ibin][icomp][iatt][itype]] = value;
  return;
}


void atmFitPars::setParameter(int ipar, double value){
  pars[ipar]=value; //< set 1D parameter array
  // update systematic and histogram modification parameters as well
  // these parameters should be kept in sync
  if (ipar>=(nTotPars-nSysPars)) sysPar[ipar-nTotPars+nSysPars] = value;
  else{
    histoPar[binOfPar[ipar]][compOfPar[ipar]][attOfPar[ipar]][typeOfPar[ipar]]=value;
  }
  return;
}

atmFitPars::atmFitPars(int isamp, int ibin, int icomp, int iatt, const char* systype){
  nSamples = isamp;
  nBins    = ibin;
  nComponents    = icomp;
  nAttributes   = iatt;;
  initPars(systype);
}

atmFitPars::atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst){
  nSamples = isamp;
  nBins    = ibin;
  nComponents    = icomp;
  nAttributes   = iatt;
  nSysPars = nsyst;
  for (int isys=0;isys<nSysPars;isys++){
    sysPar[isys]=1.;
  }
  initPars();
}


//initialize pars to preset values, the number of bins, components, attributes and systematic parameters must be previously set
void atmFitPars::initPars(const char* systype){

  //initialize histogram pars
  int index = 0;
  for (int ibin=0;ibin<nBins;ibin++){
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        histoPar[ibin][icomp][iatt][0]=1.0;
        pars[index]=1.0;
        parUnc[index]=0.005; //rough estimate of uncertainty
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=0;
        fixPar[index]=0;
        parIndex[ibin][icomp][iatt][0]=index;  
        index++;
        histoPar[ibin][icomp][iatt][1]=0.0;
        parIndex[ibin][icomp][iatt][1]=index;
        pars[index]=0.0;
        parUnc[index]=1.0;  //rough estimate of uncertainty
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=1;
        fixPar[index]=0;
        index++;
      }
    }
  }

  /////////////////////////////////////////////////
  // initialize systematic error parameters
  TString stype = systype;
  nSysPars=0;

  if (!stype.CompareTo("tn186")){
    //CCQE xsec norm bin 1//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm  bin 2//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.25;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 3//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.1;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 4//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.05;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //SubGeV flux norm//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.25;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //MultiGeV flux norm//
    sysPar[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.15;
    nSysPars++;
    //CCnQE xsec norm//
    sysPar[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.2;
    nSysPars++;
    //NC xsec norm
    sysPar[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.2;
    nSysPars++;
    //mu/e xsec ratio
    sysPar[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.05;
    nSysPars++;
  }

  // no parameters
  if (!stype.CompareTo("none")){
    nSysPars=0;
  }

  // cosmic muon parameters
  if (!stype.CompareTo("cosmic")){
    //FV Bin 0 norm
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.10;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //FV Bin 1 norm
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.10;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //FV Bin 2 norm
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.10;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
  }


  //add systematics to 1D parameter arrays
  for (int isys=0;isys<nSysPars;isys++){
    pars[index]=sysPar[isys];
    parUnc[index]=sysParUnc[isys];
    sysParIndex[isys] = index;
    index++;
  }

  // initialize normalization parameters
  int normindex = index;
  int normpars  = 0;
  for (int ibin=0; ibin<nBins; ibin++){
    for (int isamp=0; isamp<nSamples; isamp++){
       histoNorm[isamp][ibin] = 1.; //< default histo norm is one
       pars[normindex] = 1.0;
       parUnc[normindex] = 0.1;
       sysParUnc[nSysPars] = 0.1;
       sysParDefault[nSysPars] = 1.0;
       normParIndex[isamp][ibin] = normindex;
       normindex++;
       //count these as systematic parameters if using in fit
       if (flgUseNormPars){
         nSysPars++;
	 index++;
       }
    }
  }

  // end systemaitc parmeter initializations
  /////////////////////////////////////////////////


  ////////////////////////////////////////////
  // fix total number of parameters and print values
  nTotPars = index;
  cout<<"Total number of fit parameters: "<<nTotPars<<endl;
  for (int kpar=0;kpar<nTotPars;kpar++){
    cout<<"par "<<kpar<<" value: "<<pars[kpar]<<endl;
  }
  return;
}



          
#endif
