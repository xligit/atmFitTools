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
  atmFitPars(const char* parfile); //constructs from parameter file

  ///////////////////////////////////////////////////////////////
  //numbers of various parametrs
  int nSamples;
  int nBins;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nTotPars;

  ///////////////////////////////////////////////////////////////////
  //parameter values
  double histoPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncLo[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  double sysPar[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  double pars[4000];
  double parUnc[4000];
  int   fixPar[4000]; //< array of fix flags for parameters
  double bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2];
  double norm;  

  //////////////////////////////////////////////////////////////
  //methods
  void setNorm(double x){norm=x;}
  void initPars(const char* systype=""); //< sets parameters to initial values
  int getParIndex(int ibin, int icomp, int iatt, int imod){return parIndex[ibin][icomp][iatt][imod];}
  double getParameter(int ipar){return pars[ipar];}
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
//set parameters back to defaults
void atmFitPars::resetDefaults(){

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

  //initialize systematic error parameters
  for (int isyst=0;isyst<nSysPars;isyst++){
    pars[index]=1.0;
    index++;
    sysPar[isyst]=1.0;
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

//construct from parameter file
atmFitPars::atmFitPars(const char* parfilename){
  /////////////////////////////////////
  //fill shared parameters from file
  cout<<"atmFitPars: reading parameter file: "<<parfilename<<endl;
  sharedPars sharedpar(parfilename);
  nSamples = sharedpar.nSamples;
  cout<<"  nSamples: "<<nSamples<<endl;
  nComponents = sharedpar.nComponents;
  cout<<"  nComponents: "<<nComponents<<endl;
  nBins = sharedpar.nFVBins;
  cout<<"  nBins: "<<nBins<<endl;
  nSysPars = sharedpar.nSysPars;
  cout<<"  nSysPars: "<<nSysPars<<endl;
  nAttributes = sharedpar.nAttributes;
  cout<<"  nAttributes: "<<nAttributes<<endl;

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
  pars[ipar]=value;
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

  //initialize systematic error parameters
  TString stype = systype;
  nSysPars=0;

  if (!stype.CompareTo("tn186")){
    //CCQE xsec norm bin 1//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm  bin 2//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.25;
    nSysPars++;
    //CCQE xsec norm bin 3//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.1;
    nSysPars++;
    //CCQE xsec norm bin 4//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.05;
    nSysPars++;
    //SubGeV flux norm//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.25;
    nSysPars++;
    //MultiGeV flux norm//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.15;
    nSysPars++;
    //CCnQE xsec norm//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.2;
    nSysPars++;
    //NC xsec norm
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.2;
    nSysPars++;
    //mu/e xsec ratio
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.05;
    nSysPars++;
  }

  if (!stype.CompareTo("none")){
    nSysPars=0;
  }

  //fix 1D parameter arrays
  for (int isys=0;isys<nSysPars;isys++){
    pars[index]=sysPar[isys];
    parUnc[index]=sysParUnc[isys];
    index++;
  }
  nTotPars = index;
  cout<<"Total number of fit parameters: "<<nTotPars<<endl;
  for (int kpar=0;kpar<nTotPars;kpar++){
    cout<<"par "<<kpar<<" value: "<<pars[kpar]<<endl;
  }
  return;
}



          
#endif
