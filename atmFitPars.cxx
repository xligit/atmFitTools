#ifndef ATMFITPARS_C
#define ATMFITPARS_C

#include "atmFitPars.h"

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

/////////////////////////////////////////////////
int atmFitPars::getParBin(int ipar){

  int value = -1;

  //get more information for effective parameters
  for (int ibin=0; ibin<nBins; ibin++){
    for (int icomp=0; icomp<nComponents; icomp++){
      for (int iatt=0; iatt<nAttributes; iatt++){
        for (int imod=0; imod<2; imod++){
          int parindex = getParIndex(ibin,icomp,iatt,imod);
          if (parindex==ipar){
            value = ibin;
            return value;
          }
        }
      }
    }
  }

  return value;

}

/////////////////////////////////////////////////
int atmFitPars::getParComp(int ipar){

  int value = -1;

  //get more information for effective parameters
  for (int ibin=0; ibin<nBins; ibin++){
    for (int icomp=0; icomp<nComponents; icomp++){
      for (int iatt=0; iatt<nAttributes; iatt++){
        for (int imod=0; imod<2; imod++){
          int parindex = getParIndex(ibin,icomp,iatt,imod);
          if (parindex==ipar){
            value = icomp;
            return value;
          }
        }
      }
    }
  }

  return value;

}

/////////////////////////////////////////////////
int atmFitPars::getParAtt(int ipar){

  int value = -1;

  //get more information for effective parameters
  for (int ibin=0; ibin<nBins; ibin++){
    for (int icomp=0; icomp<nComponents; icomp++){
      for (int iatt=0; iatt<nAttributes; iatt++){
        for (int imod=0; imod<2; imod++){
          int parindex = getParIndex(ibin,icomp,iatt,imod);
          if (parindex==ipar){
            value = iatt;
            return value;
          }
        }
      }
    }
  }

  return value;

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

#ifdef T2K
void atmFitPars::proposeStep()
{
  for (int i = 0; i < nTotPars - nSysPars; ++i) {
    if (!fixPar[i]) parsProp[i] = rnd->Gaus(pars[i], parUnc[i]*fScale); // random walk
  }
  for (int i = 0; i < 2; ++i) {
    if (!fixPar[i+nTotPars-nSysPars]) parsProp[i+nTotPars-nSysPars] = rnd->Gaus(sysParNom[i], sysParUnc[i]*fScale);
    while (pars[i+nTotPars-nSysPars]<0) parsProp[i+nTotPars-nSysPars] = rnd->Gaus(sysParNom[i], sysParUnc[i]*fScale);
  }
  if (!fixPar[3+nTotPars-nSysPars]) {
    if (rnd->Uniform(0,2)>1) parsProp[3+nTotPars-nSysPars] = 0;
    else parsProp[3+nTotPars-nSysPars] = 1;
  }
  cov->proposeStep();
  for (int i = 3; i < nSysPars; ++i) {
    if (!fixPar[i+nTotPars-nSysPars]) parsProp[i+nTotPars-nSysPars] = cov->getProposed(i-3);
  }
}

void atmFitPars::acceptStep()
{
  cov->acceptStep();
  for (int i = 0; i < nTotPars; ++i) pars[i] = parsProp[i];
}

atmFitPars::atmFitPars(const std::string parfilename, covBase *covm){

  /////////////////////////////////////
  //fill shared parameters from file
  cout<<"atmFitPars: reading parameter file: "<<parfilename<<endl;
  runpars = new sharedPars(parfilename.c_str());
  runpars->readParsFromFile();
  nSamples = runpars->nSamples;
  nModes = NMODE;
  cout<<"  nSamples: "<<nSamples<<endl;
  nComponents = runpars->nComponents;
  cout<<"  nComponents: "<<nComponents<<endl;
  nBins = runpars->nFVBins;
  cout<<"  nBins: "<<nBins<<endl;
  nAttributes = runpars->nAttributes;
  cout<<"  nAttributes: "<<nAttributes<<endl;
  normFactor = runpars->normFactor;
  cout<<"MC normalization: "<<normFactor<<endl;
  rnd = new TRandom3();
  fScale = 1;
  TString systype = runpars->sysParType;
  initPars(systype.Data());
  if (covm) setCov(covm);
}
#endif

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
 // nSysPars = runpars->nSysPars;
 // cout<<"atmFitPars:  nSysPars: "<<nSysPars<<endl;
  nAttributes = runpars->nAttributes;
  cout<<"atmFitPars:  nAttributes: "<<nAttributes<<endl;
  normFactor = runpars->normFactor;
  cout<<"MC normalization: "<<normFactor<<endl;
  flgUseNormPars = runpars->flgUseNormPars;
  if (flgUseNormPars) cout<<"atmFitPars: Using normalization pars"<<endl;
  TString sysType = runpars->sysParType;
  int flgFixAllSmearPars = runpars->flgFixAllSmearPars;
  if (flgFixAllSmearPars){
    fixAllSmearPars(1);
  }
  //fill all initial parameter values and count number of pars
  initPars(sysType.Data());
}



void atmFitPars::printPars(int ipar){
  if (ipar>=0){
    cout<<"PAR: "<<ipar<<" = "<<pars[ipar]<<" +/- "<<parUnc[ipar]<<endl;
    return;
  }
  else
  {
  cout<<"$$$ CURRENT PARAMETER VALUES My ASS$$$"<<endl;
  for (int i=0;i<nTotPars;i++){
    cout<<"PAR: "<<i<<" = "<<pars[i]<<" +/- "<<parUnc[i]<<endl;
  }
  cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  return;
  }
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
#ifndef T2K
  for (int i=0;i<nSysPars;i++){
     parval = randy2->Gaus(sysPar[i],(sysParUnc[i]/2.));
     if (parval<0) parval=0.;
     cout<<"par "<<i<<" is "<<parval<<endl;
     setSysParameter(i,parval);
  }
#else
  if (sysType=="t2k" || sysType=="banff") {
    parval = rnd->Gaus(sysPar[0],sysParUnc[0]/2.);
    if (parval < 0) parval = 0;
    setSysParameter(0, parval);
    parval = rnd->Gaus(sysPar[1],sysParUnc[1]/2.);
    if (parval < 0) parval = 0;
    setSysParameter(1, parval);
    parval = rnd->Uniform(0,1);
    if (parval > 0.5) setSysParameter(2, 0);
    else setSysParameter(2, 1);
    cov->proposeStep();
  }
#endif
}

/* use printPars() instead!!;;
void atmFitPars::printParValues(){
  cout<<"-------Parameter Values-------"<<endl;
  for (int ipar=0;ipar<nTotPars;ipar++){
    cout<<"par "<<ipar<<": "<<pars[ipar]<<endl;
  }
}
*/

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
  if (ipar>=(nTotPars-nSysPars)) {
    sysPar[ipar-nTotPars+nSysPars] = value;
#ifdef T2K
    if (ipar-nTotPars+nSysPars>2) cov->setPar(ipar-nTotPars+nSysPars-3, value);
#endif
  }
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
#ifdef T2K
  rnd = new TRandom3();
  fScale = 1;
#endif
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
#ifdef T2K
  rnd = new TRandom3();
  fScale = 1;
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//initialize pars to preset values, the number of bins, components, attributes and systematic parameters must be previously set
void atmFitPars::initPars(const char* systype){

  //for parameter names
  TString basename = "parameter_";


  //initialize histogram parameterss
  int index = 0;
  for (int ibin=0;ibin<nBins;ibin++){
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        // "smear" parameters
        histoPar[ibin][icomp][iatt][0]=1.0; 
        pars[index]=1.0;
        parDefaultValue[index] = 1.0;
        parUnc[index]=0.005; //rough estimate of uncertainty
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=0;
        fixPar[index]=0;
        parIndex[ibin][icomp][iatt][0]=index;  
        TString parname = basename.Data();
        parname.Append(Form("bin%d_",ibin));
        parname.Append(Form("comp%d_",icomp));
        parname.Append(Form("att%d_",iatt));
        parname.Append("smear");
        parName[index] = parname.Data();
        index++;
        // "bias" parameters
        histoPar[ibin][icomp][iatt][1]=0.0;
        parIndex[ibin][icomp][iatt][1]=index;
        pars[index]=0.0;
        parDefaultValue[index] = 0.0;
        parUnc[index]=1.0;  //rough estimate of uncertainty
        binOfPar[index]=ibin;
        compOfPar[index]=icomp;
        attOfPar[index]=iatt;
        typeOfPar[index]=1;
        parname = basename.Data();
        parname.Append(Form("bin%d_",ibin));
        parname.Append(Form("comp%d_",icomp));
        parname.Append(Form("att%d_",iatt));
        parname.Append("bias");
        parName[index] = parname.Data();
        fixPar[index]=0;
        index++;
      }
    }
  }

  /////////////////////////////////////////////////
  // initialize systematic error parameters
  TString stype = systype;
  nSysPars=0;
  nNormPars=0;
  
  if (!stype.CompareTo("tn186simple")){
    //CCQE xsec norm bin 1//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm  bin 2//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.411;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 3//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.216;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 4//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.155;
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

  else if (!stype.CompareTo("tn186")){
    //CCQE xsec norm bin 1//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 1.0;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm  bin 2//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.411;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 3//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.216;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 4//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.155;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 5//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.125;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 6//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.105;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 7//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.0805;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 8//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.066;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 9//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.0542;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 10//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.0398;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 11//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.0344;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 12//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.0226;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 13//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.0165;
    sysParDefault[nSysPars] = 1.0;
    nSysPars++;
    //CCQE xsec norm bin 14//
    sysPar[nSysPars] = 1.0;
    sysParUnc[nSysPars] = 0.00903;
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

  // debugging parameters
  if (!stype.CompareTo("debug")){
   nSysPars = 1;
   sysPar[0] = 1.0;
   sysParDefault[0] = 1.0;
   sysParUnc[0] = 0.05;
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

  if (!stype.CompareTo("t2k") || !stype.CompareTo("banff")) {
    nSysPars = 3; // two flux errors + hadron multiplicity
    // flux
    sysParNom[0] = 1.0; sysPar[0] = 1.0;    sysParUnc[0] = 0.25; // sub-GeV flux norm
    sysParNom[1] = 1.0; sysPar[1] = 1.0;    sysParUnc[1] = 0.15; // multi-GeV flux norm
    sysParUp[0] = 9999.; sysParLow[0] = 0.;
    sysParUp[1] = 9999.; sysParLow[1] = 0.;
    sysParDefault[0] = 1.0; sysParDefault[1] = 1.0;
    // hadron multiplicity
    sysParNom[2] = 0; sysPar[2] = 0; sysParDefault[2] = 0;
    sysParName[0] = "FLUX_SUB";
    sysParName[1] = "FLUX_MUL";
    sysParName[2] = "HAD_MULT";
    // xsec errors
    std::cout<<"Please set covariance matrix by calling atmFitPars::setCov(covBase *)\n"
	     <<"stype corresponds to "<<stype.Data()<<std::endl;
  }

  //add systematics to 1D parameter arrays
  for (int isys=0;isys<nSysPars;isys++){
    pars[index]=sysPar[isys];
    parDefaultValue[index] = sysParDefault[isys];
    parUnc[index]=sysParUnc[isys];
    sysParIndex[isys] = index;;
    TString parname = basename.Data();
    parname.Append(Form("_syspar%d",isys));
    parName[index] = parname.Data();
    index++;
  }

#ifndef T2K
  // initialize normalization parameters
  // and add to 1D array
  int normindex = index;
  int normpars  = 0;
  for (int ibin=0; ibin<nBins; ibin++){
    for (int isamp=0; isamp<nSamples; isamp++){
       histoNorm[isamp][ibin] = 1.; //< default histo norm is one
       pars[normindex] = 1.0;
       TString parname = basename.Data();
       parname.Append(Form("_normpar%d",nNormPars));
       parName[normindex] = parname.Data();
       parDefaultValue[normindex] = 1.0;
       parUnc[normindex] = 0.1;
       sysParUnc[nSysPars] = 0.1;
       sysParDefault[nSysPars] = 1.0;
       normParIndex[isamp][ibin] = normindex;
       normindex++;
       //count these as systematic parameters if using in fit
       if (flgUseNormPars){
         nSysPars++;
         nNormPars++;
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
#endif
  return;
}

#ifdef T2K
void atmFitPars::setCov(covBase *covariance)
{
  std::cout<<"setting covariance matrix"<<std::endl;
  cov = covariance;
  int index = nTotPars + 3;
  nSysPars += cov->getNPar();
  for (int i = 3; i < nSysPars; ++i) {
    //std::cout<<i<<std::endl;
    sysPar[i] = cov->getNominal(i-3);
    sysParNom[i] = cov->getNominal(i-3);
    sysParUnc[i] = cov->getUncertainty(i-3);
    sysParUp[i] = cov->getUp(i-3);
    sysParLow[i] = cov->getLow(i-3);
    sysParName[i] = cov->getParName(i-3);
    pars[nTotPars+i-3] = sysPar[i];
    parUnc[nTotPars+i-3] = sysParUnc[i];
    parDefaultValue[nTotPars+i-3] = cov->getInit(i-3);
    sysParIndex[i] = index;
    index++;
  }
  nTotPars = index;

  // initialize normalization parameters
  // and add to 1D array
  // this should be changed when performing a T2K-SK joint analysis
  // to account for the higher energy portion of SK events 
  // whose xsec systematics are not properly constrained by T2K ND280 measurements
  // i.e. different energy range should have different systematic parameter
  int normindex = index;
  int normpars  = 0;
  for (int ibin=0; ibin<nBins; ibin++){
    for (int isamp=0; isamp<nSamples; isamp++){
       histoNorm[isamp][ibin] = 1.; //< default histo norm is one
       pars[normindex] = 1.0;
       parDefaultValue[normindex] = 1.0;
       parUnc[normindex] = 0.1;
       sysParUnc[nSysPars] = 0.1;
       sysParDefault[nSysPars] = 1.0;
       normParIndex[isamp][ibin] = normindex;
       normindex++;
       //count these as systematic parameters if using in fit
       if (flgUseNormPars){
         nSysPars++;
         nNormPars++;
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

  cout<<"    - Total number of fit parameters now is: "<<nTotPars<<endl;
  cout<<"    - Total number of systematic parameters: "<<nSysPars<<endl;

}
#endif
          
#endif
