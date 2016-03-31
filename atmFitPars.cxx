#include "atmFitPars.h"

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

void atmFitPars::proposeStep()
{
  for (int i = 0; i < nTotPars - nSysPars; ++i) {
    if (!fixPar[i]) parsProp[i] = rnd->Gaus(pars[i], parUnc[i]*fScale); // random walk
  }
  for (int i = 0; i < 2; ++i) {
    if (!fixPar[i+nSysPars-nSysPars]) parsProp[i+nSysPars-nSysPars] = rnd->Gaus(sysParNom[i], sysParUnc[i]*fScale);
    while (pars[i+nSysPars-nSysPars]<0) parsProp[i+nSysPars-nSysPars] = rnd->Gaus(sysParNom[i], sysParUnc[i]*fScale);
  }
  if (!fixPar[3+nSysPars-nSysPars]) {
    if (rnd->Uniform(0,2)>1) parsProp[3+nSysPars-nSysPars] = 0;
    else parsProp[3+nSysPars-nSysPars] = 1;
  }
  cov->proposeStep();
  for (int i = 3; i < nSysPars; ++i) {
    if (!fixPar[i+nSysPars-nSysPars]) parsProp[i+nSysPars-nSysPars] = cov->getProposed(i-3);
  }
}

void atmFitPars::acceptStep()
{
  cov->acceptStep();
  for (int i = 0; i < nTotPars; ++i) pars[i] = parsProp[i];
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
atmFitPars::atmFitPars(const char* parfilename, covBase *covm){
  /////////////////////////////////////
  //fill shared parameters from file
  cout<<"atmFitPars: reading parameter file: "<<parfilename<<endl;
  runpars = new sharedPars(parfilename);
  runpars->readParsFromFile();
  nSamples = runpars->nSamples;
  nModes = 9;
  cout<<"  nSamples: "<<nSamples<<endl;
  nComponents = runpars->nComponents;
  cout<<"  nComponents: "<<nComponents<<endl;
  nBins = runpars->nFVBins;
  cout<<"  nBins: "<<nBins<<endl;
  //nSysPars = runpars->nSysPars;
  //cout<<"  nSysPars: "<<nSysPars<<endl;
  nAttributes = runpars->nAttributes;
  cout<<"  nAttributes: "<<nAttributes<<endl;
  rnd = new TRandom3();
  fScale = 1;
  TString systype = runpars->sysParType;
  initPars(systype.Data());
  if (covm) setCov(covm);
}

void atmFitPars::printPars(){
  cout<<"$$$ CURRENT PARAMETER VALUES $$$"<<endl;
  for (int i=0;i<nTotPars;i++){
    cout<<"PAR: "<<i<<" = "<<pars[i]<<" +/- "<<parUnc[i]<<endl;
  }
  //cov->PrintNominal();
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
}

TRandom2* randy2 = new TRandom2();

void atmFitPars::setRandSysPar(){
  double parval;
  if (sysType=="tn186") {
    for (int i=0;i<nSysPars;i++){
      parval = randy2->Gaus(sysPar[i],(sysParUnc[i]/2.));
      if (parval<0) parval=0.;
      cout<<"par "<<i<<" is "<<parval<<endl;
      setSysParameter(i,parval);
    }
  } else if (sysType=="t2k" || sysType=="banff") {
    parval = randy2->Gaus(sysPar[0],sysParUnc[0]/2.);
    if (parval < 0) parval = 0;
    setSysParameter(0, parval);
    parval = randy2->Gaus(sysPar[1],sysParUnc[1]/2.);
    if (parval < 0) parval = 0;
    setSysParameter(1, parval);
    parval = randy2->Uniform(0,1);
    if (parval > 0.5) setSysParameter(2, 0);
    else setSysParameter(2, 1);
    cov->proposeStep();
  }
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
}

void atmFitPars::fixParameter(int ipar){
  fixPar[ipar] = 1;
}

void atmFitPars::setSysParameter(int ipar, double value){
  sysPar[ipar]=value;
  pars[nTotPars-nSysPars+ipar]=value;
}

void atmFitPars::setParameter(int ibin, int icomp, int iatt, int itype, double value){
  histoPar[ibin][icomp][iatt][itype] = value;
  pars[parIndex[ibin][icomp][iatt][itype]] = value;
}

void atmFitPars::setParameter(int ipar, double value){
  pars[ipar]=value;
  if (ipar>=(nTotPars-nSysPars)) sysPar[ipar-nTotPars+nSysPars] = value;
  else{
    histoPar[binOfPar[ipar]][compOfPar[ipar]][attOfPar[ipar]][typeOfPar[ipar]]=value;
  }
}

atmFitPars::atmFitPars(int isamp, int ibin, int icomp, int iatt, const char* systype){
  nSamples = isamp;
  nBins    = ibin;
  nComponents    = icomp;
  nAttributes   = iatt;;
  initPars(systype);
  rnd = new TRandom3();
  fScale = 1;
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
  rnd = new TRandom3();
  fScale = 1;
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
  } else if (!stype.CompareTo("t2k") || !stype.CompareTo("banff")) {
    nSysPars = 3; // two flux errors + hadron multiplicity
    // flux
    sysParNom[0] = 1.0; sysPar[0] = 1.0;    sysParUnc[0] = 0.25; // sub-GeV flux norm
    sysParNom[1] = 1.0; sysPar[1] = 1.0;    sysParUnc[1] = 0.15; // multi-GeV flux norm
    sysParUp[0] = 9999.; sysParLow[0] = 0.;
    sysParUp[1] = 9999.; sysParLow[1] = 0.;
    // hadron multiplicity
    sysParNom[2] = 0; sysPar[2] = 0;
    sysParName[0] = "FLUX_SUB";
    sysParName[1] = "FLUX_MUL";
    sysParName[2] = "HAD_MULT";
    // xsec errors
    std::cout<<"Please set covariance matrix by calling atmFitPars::setCov(covBase *)\n"
	     <<"stype corresponds to "<<stype.Data()<<std::endl;
  }
  
  if (!stype.CompareTo("none")){
    nSysPars=0;
  }

  //fix 1D parameter arrays
  for (int isys=0;isys<nSysPars;isys++){
    pars[index]=sysPar[isys];
    parUnc[index]=sysParUnc[isys];
    //std::cout<<isys<<" "<<index<<std::endl;
    index++;
  }
  nTotPars = index;
  cout<<"Total number of fit parameters: "<<nTotPars<<endl;
  /*
  for (int kpar=0;kpar<nTotPars;kpar++){
    cout<<"par "<<kpar<<" value: "<<pars[kpar]<<endl;
    }*/

}

void atmFitPars::setCov(covBase *covariance)
{
  std::cout<<"setting covariance matrix"<<std::endl;
  cov = covariance;
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
  }
  nTotPars += cov->getNPar();
  cout<<"Total number of fit parameters now is: "<<nTotPars<<endl;
  cout<<"Total number of systematic parameters: "<<nSysPars<<endl;
  //cov->PrintNominal();
}
