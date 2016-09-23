#ifndef HISTOCOMPARE_C
#define HISTOCOMPARE_C

#include "histoCompare.h"
#include "markovTools.cxx"
#include <time.h>

histoCompare* histoCompare::staticthis;


/////////////////////////////////////////////
//get a rough estimation (fron 1D likelihood proflie)
//of the uncertainty of each fit parameter
//Requires a guess of parUnc[] in atmFitPars initialization
//This should be run before running MCMC, since these
//rough uncertainties are used to proposes MCMC steps
void histoCompare::calcRoughParErr(){
   cout<<"histoCompare: "<<"calculating rough uncertainties"<<endl;  
  //print final results
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    if (thePars->fixPar[ipar]) continue;
    errParLo[ipar]=getErrLo(ipar);
    errParHi[ipar]=getErrHi(ipar);
    thePars->parUnc[ipar]=(errParHi[ipar]-errParLo[ipar]);
    cout<<"  PAR "<<ipar<<" FIT RESULT: "<<thePars->pars[ipar]<<" +/- : "<<thePars->parUnc[ipar]<<endl;
  }
  return;
}

void histoCompare::readFitPars(const char* filename){
  thePars->readPars(filename);
  return;
}

void histoCompare::timetest(int ntry){
  clock_t t1,t2;
  double par[100];
  int parindex = 0;
  for (int ibin=0;ibin<nBin;ibin++){
    for (int icomp=0;icomp<nComp;icomp++){
      for (int iatt=0;iatt<nAtt;iatt++){
        for (int imod=0;imod<2;imod++){
           par[parindex]=Par[ibin][icomp][iatt][imod];
           parindex++;
        }
      }
    }
  }
  double diff;
  int itry=0;
  t1=clock();
  while (itry<ntry){
    //getTotLnL1D(result,par); 
    getTotLnL();
    itry++;
  }
  t2=clock();
  diff = ((double)t2-(double)t1)/(double)ntry;
  double diff1 = diff;
  cout<<"time 1: "<<diff<<endl;
  itry=0;
  t1=clock();
  while (itry<ntry){
    //getTotLnL1D(result,par); 
    getTotLnL();
    itry++;
  }
  t2=clock();
  diff = ((double)t2-(double)t1)/(double)ntry;
  cout<<"time 2: "<<diff<<endl;
  cout<<"changes: "<<diff1-diff<<endl;
  return;
}


void histoCompare::tuneDEMCMC(int ncycles,int nsteps, double goal){
 
  double result = 0.;
  markovTools* mcmc = new markovTools(thePars); //< create markovTools object
  mcmc->setTuneParameter(tunePar);
  mcmc->setDiffChain(diffChainFileName.Data());

  //set initial state
  result = getTotLnL();
  mcmc->setL(result);//< sets the initial likelihood
  double Linit = result;

  //run tuning
  for (int icycle=0;icycle<ncycles;icycle++){
    double xaccepted=0.;
    int   istep=0;
    while (istep<nsteps){
      cout<<"----------"<<"step : "<<istep<<"---------------"<<endl;;
      mcmc->proposeDifferentialStep(); //< propose a new step using differential mcmc
      result = getTotLnL();
      cout<<"hc: Likelihood "<<Linit<<" -> "<<result<<" diff: "<<result-Linit<<endl;
      if (mcmc->acceptStepLnL(result)){ //< check if new step is accepted
        xaccepted++; 
      }
      istep = mcmc->iStep;
    }

    double rate = xaccepted/(double)nsteps;
    cout<<"acceptance rate: "<<rate<<endl;
    cout<<"tune parameter: "<<tunePar<<endl;
    if ((rate>22.0)&&(rate<30.0)) return;
    tunePar*=(rate/goal);
    cout<<"new tune parameter: "<<tunePar<<endl;;
  }
  return; 
}



void histoCompare::tuneMCMC(int ncycles,int nsteps,double goal){
 
  double result = 0.;
// markovTools* mc = new markovTools(npars); //< create markovTools object
  markovTools* mcmc = new markovTools(thePars); //< create markovTools object
  mcmc->setTuneParameter(tunePar);


  //set initial state
  result = getTotLnL();
  mcmc->setL(result);//< sets the initial likelihood
  double Linit = result;

  //run tuning
  for (int icycle=0;icycle<ncycles;icycle++){
    double xaccepted=0.;
    int   istep=0;
    while (istep<nsteps){
      cout<<"----------"<<"step : "<<istep<<"---------------"<<endl;;
      mcmc->proposeStep(); //< propose a new step
      result = getTotLnL();
      cout<<"hc: Likelihood "<<Linit<<" -> "<<result<<" diff: "<<result-Linit<<endl;
      if (mcmc->acceptStepLnL(result)){ //< check if new step is accepted
        xaccepted++; 
      }
      istep = mcmc->iStep;
    }

    double rate = xaccepted/(double)nsteps;
    cout<<"acceptance rate: "<<rate<<endl;
    cout<<"tune parameter: "<<tunePar<<endl;
    tunePar*=(rate/goal);
    cout<<"new tune parameter: "<<tunePar<<endl;;
  }
  return; 
}


void histoCompare::tuneMCMCOld(int ncycles,int nsteps,double goal){
 
//  //setup mc mcmctools
  const int npars = thePars->nTotPars; //< total number of parameters in fit
  double par[npars]; //< container for parameters
  int parindex = 0;
  double result = 0.;
 markovTools* mc = new markovTools(npars); //< create markovTools object
//  markovTools* mc = new markovTools(thePars); //< create markovTools object
  mc->setTuneParameter(tunePar);

  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    par[ipar]=thePars->getParameter(ipar); //< parameter array
    mc->setParVar(ipar,thePars->parUnc[ipar]); //< set parameter variance
    mc->setFixPar(ipar,thePars->fixPar[ipar]); //< set parameter fix flag
  }

  //set initial state
 // result = getTotLnL();
  getTotLnL1D(result,npars, par);//< get total likelihood from 1D array
  mc->setL(result);//< sets the initial likelihood
  double Linit = result;

  //run tuning
  for (int icycle=0;icycle<ncycles;icycle++){
    double xaccepted=0.;
    int   istep=0;
    while (istep<nsteps){
      istep = mc->iStep;
      mc->proposeStep(par); //< propose a new step
      getTotLnL1D(result, npars,par);  //< get likelihood of new step
      //result = getTotLnL();
      cout<<result-Linit<<endl;
      if (mc->acceptStepLnL(result,par)){ //< check if new step is accepted
        xaccepted++; 
      }
    }

    double rate = xaccepted/(double)nsteps;
    cout<<"acceptance rate: "<<rate<<endl;
    cout<<"tune parameter: "<<tunePar<<endl;
    tunePar*=(rate/goal);
    cout<<"new tune parameter: "<<tunePar<<endl;
  }
  return; 
}

///////////////////////////////////////////////
//Run a MCMC of length nsteps
void histoCompare::runDiffMCMC(int nsteps){

  ///////////////////////////////////////
  //setup mcmc tools
  markovTools* mc = new markovTools(thePars);
  mc->setTuneParameter(tunePar);

  ///////////////////////////////////////////////
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    mc->setParVar(ipar,thePars->parUnc[ipar]);
  }
  
  ///////////////////////////////////////////////////
  //set initial state
  double result = getTotLnL();
  mc->setL(result);//< sets the initial likelihood

  //loop through steps
  int currentstep=0;
  while (currentstep<nsteps){
    currentstep = mc->iStep;
    //mc->proposeStep(par); //< fills par[] with proposed params
    mc->proposeStep();
    //getTotLnL1D(result, npars,par);   
    result = getTotLnL();
    //mc->acceptStepLnL(result,par); //< if step is accepted, istep++, and written to histo
    mc->acceptStepLnLDiff(result);
  }

  ///////////////////////
  //save results
  mc->savePath();

  //////////////////////////
  return;
}


///////////////////////////////////////////////
//Run a MCMC of length nsteps
void histoCompare::runDEMCMC(int nsteps){

  ///////////////////////////////////////
  //setup mcmc tools
  markovTools* mc = new markovTools(thePars);
  mc->setTuneParameter(tunePar);
  mc->setDiffChain(diffChainFileName.Data());

  ///////////////////////////////////////////////
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    mc->setParVar(ipar,thePars->parUnc[ipar]);
  }
  
  ///////////////////////////////////////////////////
  //set initial state
  double result = getTotLnL();
  mc->setL(result);//< sets the initial likelihood

  //loop through steps
  int currentstep=0;
  while (currentstep<nsteps){
    currentstep = mc->iStep;
    mc->proposeDifferentialStep();
    result = getTotLnL();
    mc->acceptStepLnL(result);
  }

  ///////////////////////
  //save results
  mc->savePath();

  //////////////////////////
  return;
}




///////////////////////////////////////////////
//Run a MCMC of length nsteps
void histoCompare::runMCMC(int nsteps){

  ///////////////////////////////////////
  //setup mcmc tools

  markovTools* mc = new markovTools(thePars);
  mc->setTuneParameter(tunePar);

  ///////////////////////////////////////////////
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    mc->setParVar(ipar,thePars->parUnc[ipar]);
  }
  
  ///////////////////////////////////////////////////
  //set initial state
  double result = getTotLnL();
  mc->setL(result);//< sets the initial likelihood

  //loop through steps
  int currentstep=0;
  while (currentstep<nsteps){
    currentstep = mc->iStep;
    //mc->proposeStep(par); //< fills par[] with proposed params
    mc->proposeStep();
    //getTotLnL1D(result, npars,par);   
    result = getTotLnL();
    //mc->acceptStepLnL(result,par); //< if step is accepted, istep++, and written to histo
    mc->acceptStepLnL(result);
  }

  ///////////////////////
  //save results
  mc->savePath();

  //////////////////////////
  return;
}



void histoCompare::saveFitPars(const char* filename){
  thePars->savePars(filename);
  return;
}

double histoCompare::getErrHi(int ipar){
  //scan through log likelihood by decreasing parameter until threshold is reached
  double Lthresh = 1.0; //likelihood threshold
  double Ldiff = 0.;  //difference in likelihood from current value
  double Lbest = getTotLnL(); //current likelihood value
 // cout<<"lbest: "<<Lbest<<endl;
  double parbest = thePars->pars[ipar];
  double parval = thePars->pars[ipar]; 
  double dpar = thePars->parUnc[ipar]/10.;
  double hierr;
  int ntry = 0;
  int ntrymax=5000;
  //coarse search
  while (Ldiff<Lthresh){
  //  cout<<"oldpar: "<<parval<<endl;
    parval+=dpar;
 //   cout<<"newpar: "<<parval<<endl;
    thePars->setParameter(ipar,parval); //modify parameter
 //   double Lnew = getTotLnL();
 //   cout<<"Lnew: "<<Lnew<<endl;
    Ldiff = TMath::Abs(Lbest-getTotLnL()); //check L difference;
 //   cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
    if (ntry>ntrymax) break;
  }
//  cout<<"ntry: "<<ntry<<endl;
  parval-=dpar;
  thePars->setParameter(ipar,parval);
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  //fine search
  while ((Ldiff<1)){
    parval+=dpar;
    thePars->setParameter(ipar,parval); //modify paramete
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
    ntry++;
    if (ntry>ntrymax) break;
  }
  parval-=dpar;
  thePars->setParameter(ipar,parval);
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  //very fine search
//  while ((Ldiff<1)&&(ntry<ntrymax)){
//    parval+=dpar;
//    thePars->setParameter(ipar,parval); //modify paramete
//    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
//    ntry++;
//  }
//  cout<<"ntry: "<<ntry<<endl;
  hierr = thePars->pars[ipar];
  thePars->setParameter(ipar,parbest);
  return hierr; 
}



double histoCompare::getErrLo(int ipar){
  //scan through log likelihood by decreasing parameter until threshold is reached
  double Lthresh = 1.0; //likelihood threshold
  double Ldiff = 0.;  //difference in likelihood from current value
  double Lbest = getTotLnL(); //current likelihood value
  double parbest = thePars->pars[ipar]; //current parameter value
  double parval = thePars->pars[ipar]; //floating value for estimation 
  double dpar = thePars->parUnc[ipar]/10.; //how much parameter should change between steps
//  cout<<"dpar: "<<dpar<<endl;
  double loerr;
  int ntry = 0;
  int ntrymax=5000;
  //coarse search
  while (Ldiff<Lthresh){
    parval-=dpar;
    thePars->setParameter(ipar,parval); //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
 //   cout<<"ntry: "<<ntry<<endl;
 //   cout<<"Ldiff: "<<Ldiff<<endl;
 //   cout<<"par: "<<thePars->pars[ipar]<<endl;
    ntry++;
    if (ntry>ntrymax) break;
  }
//  cout<<"ntry: "<<ntry<<endl;
  parval+=dpar;
  thePars->setParameter(ipar,parval);
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  //fine search
  while ((Ldiff<1)){
    parval-=dpar;
    thePars->setParameter(ipar,parval); //modify paramete
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
    ntry++;
    if (ntry>ntrymax) break;
  }
//  cout<<"ntry: "<<ntry<<endl;
  parval+=dpar;
  thePars->setParameter(ipar,parval);
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  //very fine search
//  while ((Ldiff<1)&&(ntry<ntrymax)){
//    parval-=dpar;
//    thePars->setParameter(ipar,parval); //modify paramete
//    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
//    ntry++;
//    if (ntry>ntrymax) break;
//  }
//  cout<<"ntry: "<<ntry<<endl;
  loerr = thePars->pars[ipar];
  thePars->setParameter(ipar,parbest);
  return loerr; 
}



void histoCompare::profileL(int ipar, double range, int npts, int sameflg){
  TString pname = "p";
  pname.Append("_profile.png");
  int iprofile =0;
  if (sameflg) iprofile = 1;
  double bestpoint = thePars->pars[ipar];
  cout<<"Fitted value: "<<bestpoint<<endl;
  double dx = range/(double)npts;
  double xx = bestpoint - (range/2.);
  double ll;
  double lbest = getTotLnL();
  if (hProf[iprofile]!=NULL) hProf[iprofile]->Delete();
  hProf[iprofile] = new TH1D("hprof","hprof",npts,xx,(xx+range));
  cout<<"histoCompare: profiling parameter "<<ipar<<endl;
  for (int ipoint=0;ipoint<npts;ipoint++){
    thePars->setParameter(ipar,xx);
    ll = getTotLnL();
    hProf[iprofile]->SetBinContent(ipoint+1,ll-lbest);
    xx+=dx;
  } 
  //set parameter back to initial value
  thePars->setParameter(ipar,bestpoint);
  hProf[iprofile]->SetLineColor(9);
  hProf[iprofile]->GetYaxis()->SetTitle("#Delta #chi^{2}");
  TString xname = "parameter ";
  xname.Append(Form("%d",ipar));
  hProf[iprofile]->GetXaxis()->SetTitle(xname.Data());
  hProf[iprofile]->Draw("c");
  if (sameflg) hProf[iprofile]->SetLineColor(kRed);
  if (sameflg){ 
    hProf[0]->Draw("c");
    hProf[iprofile]->Draw("samec");
  }
  else{
    hProf[0]->Draw("c");
  }
  cc->Print(pname.Data());
  return;
}




void histoCompare::profileL(int ibin, int icomp, int iatt, int imod, double range, int npts){
  TString pname = "p";
  pname.Append("_profile.png");
  //double bestpoint = Par[ibin][icomp][iatt][imod];
  double bestpoint = thePars->getHistoParameter(ibin,icomp,iatt,imod);


  cout<<"Fitted value: "<<bestpoint<<endl;
  double dx = range/(double)npts;
  double xx = bestpoint - (range/2.);
  double ll;
  double lbest = getTotLnL();
  if (hProf[0]) hProf[0]->Delete();
  hProf[0] = new TH1D("hprof","hprof",npts,xx,(xx+range));
  for (int ipoint=0;ipoint<npts;ipoint++){
   // cout<<"filling point" <<ipoint<<endl;
    //Par[ibin][icomp][iatt][imod] = xx;
    thePars->setParameter(ibin,icomp,iatt,imod,xx);
    ll = getTotLnL();
    hProf[0]->SetBinContent(ipoint+1,ll-lbest);
 //   hProf[0]->SetBinContent(ipoint+1,ll);

    xx+=dx;
  } 
  //set parameter back to initial value
  thePars->setParameter(ibin,icomp,iatt,imod,bestpoint);
  hProf[0]->SetLineColor(9);
  hProf[0]->Draw("c");
  cc->Print(pname.Data());
  return;
}

void histoCompare::showFitPars(int ibin,int iatt,int imod){
  const int nbinstot = nComp;
  double X[nbinstot];
  double EXL[nbinstot];
  double EXH[nbinstot];
  double EYL[nbinstot];
  double EYH[nbinstot];
  double Y[nbinstot];
  int parindex;
  if (gPar) gPar->Delete();
  hPar  = new TH1D("hpar","hpar",nbinstot,0,nbinstot);
  hParErrLo = new TH1D("hparerrlo","hparerrlo",nbinstot,0,nbinstot);
  hParErrHi = new TH1D("hparerrhi","hparerrhi",nbinstot,0,nbinstot);
  for (int icomp=0;icomp<nComp;icomp++){
    parindex=thePars->getParIndex(ibin,icomp,iatt,imod);
    hPar->SetBinContent(icomp+1,thePars->getParameter(parindex));
    X[icomp]=(double)hPar->GetBinCenter(icomp+1);
    Y[icomp]=(double)thePars->getParameter(parindex);
 //   hPar->GetXaxis()->SetBinLabel(icomp+1,parName[ibin][icomp][iatt][imod].Data());
    EYL[icomp]=(double)thePars->getParameter(parindex) -errParLo[parindex];
    EYH[icomp]=errParHi[parindex] - (double)thePars->getParameter(parindex);
    hParErrLo->SetBinContent(icomp+1,(double)EYL[icomp]);
    hParErrHi->SetBinContent(icomp+1,(double)EYH[icomp]);
    EXL[icomp]=0.4;
    EXH[icomp]=0.4;
  }
  gPar = new TGraphAsymmErrors(nbinstot,X,Y,EXL,EXH,EYL,EYH);
  gPar->SetFillColor(6);
  gPar->SetLineColor(6);
  gPar->SetMarkerStyle(8);
  gPar->Draw("a2");
  gPar->Draw("p");
  return;
}

void histoCompare::showModHiso(int isamp,int ibin, int icomp, int iatt, double smear, double bias){
  double biastmp = thePars->getHistoParameter(ibin,icomp,iatt,0);
  double smeartmp = thePars->getHistoParameter(ibin,icomp,iatt,1);
  thePars->setParameter(ibin,icomp,iatt,0,smear);
  thePars->setParameter(ibin,icomp,iatt,1,bias);
  hMod = hManager->getModHistogram(isamp,ibin,icomp,iatt); //gets the modified histogram
  hTmp = hManager->getHistogram(isamp,ibin,icomp,iatt);
//  if (hMod) hMod->Delete();
//  if (hTmp) hTmp->Delete();
//  if (hTot) hTot->Delete();
//  hMod = (TH1D*) hManager->hMC[isamp][ibin][icomp][iatt]->Clone("hclonetmp");
//  hTot = (TH1D*) hManager->hMC[isamp][ibin][icomp][iatt]->Clone("htottmp");
//  hTmp = (TH1D*) hManager->hMC[isamp][ibin][icomp][iatt]->Clone("htmp");
//  smearHisto((*hManager->hMC[isamp][ibin][icomp][iatt]),(*hMod),smear,bias);
  
//  for (int icomp=1;i<nComp;icomp++){
//     smearHisto((*hManager->hMC[isamp][ibin][icomp][iatt]),(*hTmp),smear,bias);
//     hMod->Add(hTmp);
//     hTot->Add(hManager->hMC[isamp][ibin][icomp][iatt]);
//  }
//  hMod->Rebin(rebinFactor);
 // hTot->Rebin(rebinFactor);
  hMod->SetLineColor(kRed);
  hTmp->Draw();
  hMod->Draw("sameh");
  //set parameter back to original
  thePars->setParameter(ibin,icomp,iatt,0,smeartmp);
  thePars->setParameter(ibin,icomp,iatt,1,biastmp);

  return;
}


void histoCompare::showFitDiff(int isamp,int ibin,int iatt){

  // get (normalized) histogram with parameter modifications
  hMod = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt)->Clone("hmod");

  // draw MC histograms
  hMod->SetLineColor(kBlue);
  // subtract data from modified expectation
  hMod->Add(hManager->hData[isamp][ibin][iatt],-1.);
  hMod->Draw("h");
 
  //
  return;
}



void histoCompare::showFitResult(int isamp,int ibin,int iatt){

  // get (normalized) histogram with parameter modifications
  hMod = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt)->Clone("hmod");

  // get nominal histgram (also normalized)
  hTmp = hManager->getSumHistogram(isamp,ibin,iatt);

  // draw MC histograms
  hMod->SetLineColor(kBlue);
  hTmp->SetLineColor(kRed);
 // hTmp->Draw("h");
  hMod->Draw("h");
 
  // draw data histograms
  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->Draw("samee");

  //
  return;
}

//Show the effect of varying a single set of bias and smear parameters
void histoCompare::showFitEffect(int isamp,int ibin,int icomp,int iatt){

  // get the modified version of the specified histogram
  hMod = hManager->getModHistogram(isamp,ibin,icomp,iatt); //gets the modified histogram

  // get the default sum of histograms
  hTmp = hManager->getSumHistogram(isamp,ibin,iatt); //get the sum histogram

  // add default histograms to modified histograms
  for (int jcomp=0;jcomp<nComp;jcomp++){
    if (jcomp!=icomp){
      hMod->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
    }
  }

  // modified histogram is blue
  hMod->SetLineColor(kBlue);

  // default is red
  hTmp->SetLineColor(kRed);

  // scale to data
  double thenorm = hManager->normFactor;
  hMod->Scale(thenorm);
//  hTmp->Scale(thenorm);

  // draw to same canvas
  hMod->Draw("h");
  hTmp->Draw("sameh");
  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->Draw("samee");
  return;
}

///////////////////////////////////////////////////////////////////////////////////
//Plots the specified histogram before and after modifications, as well as the data
void histoCompare::showFitHisto(int isamp,int ibin,int icomp,int iatt){
  double smear = thePars->getHistoParameter(ibin,icomp,iatt,0);
  double bias  = thePars->getHistoParameter(ibin,icomp,iatt,1);
  cout<<"SMEAR: "<<smear<<endl;
  cout<<"BIAS:  "<<bias<<endl;
  hMod = hManager->getModHistogramMC(isamp,ibin,icomp,iatt);
  hMod->SetLineColor(kBlue);
  hMod->Draw("eh");
  hTmp = hManager->getHistogram(isamp,ibin,icomp,iatt);
  hTmp->SetLineColor(kRed);
  hTmp->Draw("sameeh");
  return;
}

void histoCompare::lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg){

  for (int ipar=0;ipar<staticthis->thePars->nTotPars;ipar++){
    staticthis->thePars->setParameter(ipar,par[ipar]);
  }

  result = (double)staticthis->getTotLnL();
}


////////////////////////////////////////////
// Pre-fit parameters using MINUIT
void histoCompare::LnLPreFit(){

  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //threshold to determine if a parameter is fit or not
  //if the MC histograms for this parameter have a size less than this value,
  //don't bother fitting them!
  double nthresh = 100.;

  //sets the precision of the fits
  double parerr = 0.05;  
  
  //individually fit each parameter
  int parindex =0;

  //parameter name container
  TString parnametmp;  

  //total number of parameters to be fit
  int npars = thePars->nTotPars;

  cout<<"$$$$$$$$$$$$$$$ LNL PRE FIT $$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //fix parameters with too few events to be fit
  parindex = 0;
  for (int ibin=0;ibin<nBin;ibin++){
    for (int iatt=0;iatt<nAtt;iatt++){
      for (int icomp=0;icomp<nComp;icomp++){
        //name parameters
        parnametmp = Form("par_%d_",parindex);
        parnametmp.Append(binName[ibin].Data());
        parnametmp.Append("_");
        parnametmp.Append(compName[icomp].Data());
        parnametmp.Append("_");
        parnametmp.Append(attName[iatt].Data());
        parnametmp.Append("_");
        parnametmp.Append("smear");
        parName[ibin][icomp][iatt][0]=parnametmp.Data();
        parnametmp = Form("par_%d_",(parindex+1));
        parnametmp.Append(binName[ibin].Data());
        parnametmp.Append("_");
        parnametmp.Append(compName[icomp].Data());
        parnametmp.Append("_");
        parnametmp.Append(attName[iatt].Data());
        parnametmp.Append("_");
        parnametmp.Append("bias");
        parName[ibin][icomp][iatt][1]=parnametmp.Data();
        //get summed histogram
        hTot = (TH1D*)hManager->hMC[0][ibin][icomp][iatt]->Clone("htot");
        for (int isamp=1;isamp<nSamp;isamp++){
          hTot->Add(hManager->hMC[isamp][ibin][icomp][iatt]);
        }
        cout<<"total entries: "<<hTot->GetEntries()<<endl;
        //check to make sure histogram is above fitting threshold
        if (hTot->GetEntries()<nthresh){
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][0].Data()<<" (# ENTRIES TOO LOW!) "<<endl; 
          fixPar[ibin][icomp][iatt][0]=1;
          thePars->fixParameter(ibin,icomp,iatt,0);
          thePars->fixParameter(ibin,icomp,iatt,1);
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][1].Data()<<" (# ENTRIES TOO LOW!) "<<endl; 
        }
        hTot->Delete();
        parindex+=2;
      }
    }
  }
  cout<<"  ----------------------------------------  "<<endl;

  //setup the fitter!
  TFitter* fit = new TFitter(npars);

  //shut fitter up
  {
    double pp = -1;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }

  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //set parameters to inital values
  TString pname;
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    pname = Form("parameter%d",ipar);
    fit->SetParameter(ipar,pname.Data(),thePars->pars[ipar],parerr,0,0);
  }

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  
  parindex=0;
  //run individual bias fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
          parindex++; //< starts on odd parameter (bias only)
          //release biasparameter to be fit
          if (thePars->checkFixFlg(jbin,jcomp,jatt,1)==1){
            parindex++; //< do nothing if parameter is fixed
            continue;
          }
          fit->ReleaseParameter(parindex);
          cout<<"fitting "<<jbin<<jcomp<<jatt<<1<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,1)<<endl;
          fit->ExecuteCommand("MIGRAD",0,0);
          fit->FixParameter(parindex);
          parindex++;
      }
    }
  }
  parindex=0;

  //run individual smear fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
          //starts on even parameter (smear only)
          //release smear parameter to be fit
          if (thePars->checkFixFlg(jbin,jcomp,jatt,0)==1){
            parindex++; //<do nothing if parameter is fixed
            continue;
          }
          fit->ReleaseParameter(parindex);
          cout<<"fitting "<<jbin<<jcomp<<jatt<<0<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,0)<<endl;
          fit->ExecuteCommand("SIMPLEX",0,0);
          fit->FixParameter(parindex);
          parindex+=2;
      }
    }
  }

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }

  //release and fit flux and xsec parameters
  parindex = thePars->nTotPars-thePars->nSysPars;
  for (int isyspar=0;isyspar<thePars->nSysPars;isyspar++){
    if ((thePars->fixPar[parindex])!=1)fit->ReleaseParameter(parindex);
    parindex++;
  }
  fit->ExecuteCommand("SIMPLEX",0,0);
  parindex = 0;
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  //print final results
  for (int ipar=0;ipar<npars;ipar++){
    thePars->setParameter(ipar,fit->GetParameter(ipar));
    cout<<"  PAR "<<ipar<<" FIT RESULT: "<<thePars->pars[ipar]<<endl;
  }
  cout<<"$$$$$$$$$$$$$ END LNL PRE FIT $$$$$$$$$$$$$"<<endl;

  return;
}


void histoCompare::sysParFit(){
  //setup static this so wrapper doesn't segfault
  staticthis = this;
  //sets the precision of the fits
  double parerr = 0.001;  
  cout<<"$$$$$$$$$$$$$$$$$ LNL SINLGE PAR FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //setup the fitter!

  //total number of parameters to be fit
  int npars = thePars->nTotPars;
  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = 0;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //setup parameters
  TString aname;
  for (int ipar=0;ipar<npars;ipar++){
    aname = "parameter";
    aname.Append(Form("_%d",ipar));
    fit->SetParameter(ipar,aname.Data(),thePars->pars[ipar],parerr,0,0);
  }
 
  //fix all params
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  //release single parameter to fit
  for (int jpar=(thePars->nTotPars-thePars->nSysPars);jpar<thePars->nTotPars;jpar++){
    fit->ReleaseParameter(jpar);
  }
  //fit that thang
  fit->ExecuteCommand("SIMPLEX",0,0); 
  fit->ExecuteCommand("SIMPLEX",0,0); 
  //print results
  for (int jpar=(thePars->nTotPars-thePars->nSysPars);jpar<thePars->nTotPars;jpar++){
    cout<<"PAR: "<<jpar<<" "<<fit->GetParameter(jpar)<<endl;
  }

  return;
}


/////////////////////////////////////////////////
// Fits only parameter "ipar"
void histoCompare::singleParFit(int ipar){
  //setup static this so wrapper doesn't segfault
  staticthis = this;
  //sets the precision of the fits
  double parerr = 0.001;  
  cout<<"$$$$$$$$$$$$$$$$$ LNL SINLGE PAR FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  PARAMETER: "<<ipar<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //setup the fitter!

  //total number of parameters to be fit
  int npars = thePars->nTotPars;
  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = 0;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //setup parameters
  TString aname;
  for (int kpar=0;kpar<npars;kpar++){
    aname = "parameter";
    aname.Append(Form("_%d",kpar));
    fit->SetParameter(kpar,aname.Data(),thePars->pars[kpar],parerr,0,0);
  }
 
  //fix all params
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  //release single parameter to fit
  fit->ReleaseParameter(ipar);
  //fit that thang
  fit->ExecuteCommand("SIMPLEX",0,0); 
  fit->ExecuteCommand("SIMPLEX",0,0); 
  //print results
  cout<<"RESULT: "<<fit->GetParameter(ipar)<<endl;
 
  return;
}


void histoCompare::printFitResults(const char* directory){
  for (int isamp=0;isamp<nSamp;isamp++){
    for (int ibin=0;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
        showFitResult(isamp,ibin,iatt);
        TString plotname = directory;
        plotname.Append(hTmp->GetName());
        plotname.Append(".png");
        cc->Print(plotname.Data()); 
      }
    }
  }
  return;
}



/////////////////////////////////////////////////////////////
// Fits all parameters using MINUIT
void histoCompare::LnLFit(){
  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //sets the precision of the fits
  double parerr = 0.0001;  
  
  //individually fit each parameter
  int parindex =0;

  //parameter name container
  TString parnametmp;  

  //will we fix all the smearing parameters?
  if (flgFixAllSmearPars){
    thePars->fixAllSmearPars();
  }

  //total number of parameters to be fit
  int npars = thePars->nTotPars;

  /////////////////////////////////////////////////////////
  //run the prefit
  LnLPreFit();

  ////////////////////////////////////////////////////////
  //setup the fitter!
  cout<<"$$$$$$$$$$$$$$$$$ LNL FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = 0;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  
  //specify function to be fit
  fit->SetFCN(lnLWrapper);
  //setup parameters
  TString aname;
  for (int ipar=0;ipar<npars;ipar++){
    aname = "parameter_";
    aname.Append(Form("%d",ipar));
    fit->SetParameter(ipar,aname.Data(),thePars->pars[ipar],parerr,0,0);
  } 
  parindex = 0;

  ///////////////////////////////////////////////////////////////////
  //do individual fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){

      //start of fit block//
      //
      //fix all parameters
      for (int jpar=0;jpar<npars;jpar++){
        fit->FixParameter(jpar);
      }

      //for bin 6, fit the xsec and flux parameters as well//
      if (jbin==5){
        //release all flux and xsec pars
        for (int isyspar=(thePars->nTotPars-thePars->nSysPars);isyspar<thePars->nSysPars;isyspar++){
          if ((thePars->fixPar[parindex])!=1)fit->ReleaseParameter(isyspar);  
          //fit->ReleaseParameter(isyspar);
        }
        //fit these pars first
        fit->ExecuteCommand("SIMPLEX",0,0);
      }

      //release bias parameters
      for (int jcomp=0;jcomp<nComp;jcomp++){
        if (thePars->checkFixFlg(jbin,jcomp,jatt,1)!=1){
          fit->ReleaseParameter(thePars->getParIndex(jbin,jcomp,jatt,1));
        }
      }

      //run fit for bias parameters
      fit->ExecuteCommand("SIMPLEX",0,0);

      //now release smear parameters
      for (int jcomp=0;jcomp<nComp;jcomp++){
        if (thePars->checkFixFlg(jbin,jcomp,jatt,0)!=1){
          cout<<"fitting parameter: "<<jbin<<jcomp<<jatt<<0<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,0)<<endl;
          fit->ReleaseParameter(thePars->getParIndex(jbin,jcomp,jatt,0));
        }
      }

      fit->ExecuteCommand("SIMPLEX",0,0); //run the fit for ALL parameters

      //end of fit block//
      
    }
  }

  //set final results
  for (int ipar=0;ipar<npars;ipar++){
    thePars->setParameter(ipar,fit->GetParameter(ipar));
  }

  //calculate rough errors and print results
  calcRoughParErr();

  //end
  cout<<"  ----------------------------------------  "<<endl;
  cout<<"$$$$$$$$        FIT  COMPLLETE      $$$$$$$$"<<endl;

  return;


}

void histoCompare::getTotLnL1D(double& result,int npar, double par[]){

  for (int ipar=0;ipar<npar;ipar++){
    thePars->setParameter(ipar,par[ipar]);
  }

  result = getTotLnL();
}



////////////////////////////////////////////////
//Compute the total log liklihood by comparing all histograms
///*
double histoCompare::getTotLnL(){

  double totL = 0.;

  nDOF = 0.;
  ////////////////////////////////////////
  //contribution from histogram comparison  
//  for (int isamp=0;isamp<nSamp;isamp++){
  for (int isamp=0;isamp<1;isamp++){
    for (int ibin=5;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
       	TH1D* hPrediction = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt,1); //< get normalized histogram.
      	TH1D* hDataTmp = (TH1D*)hManager->getHistogramData(isamp,ibin,iatt);
        double partialLnL = hManager->histoLogL;
	      totL+=partialLnL;
        nDOF+=hManager->nDOF;
      }
    }
  }

  //////////////////////////////////////////////
  //contribution from flux/xsec priors
  double pull;
#ifndef T2K
  for (int isys=0;isys<thePars->nSysPars;isys++){
    pull = thePars->getSysParameter(isys)-1.;
    pull/=thePars->sysParUnc[isys];
    totL+=(0.5)*(pull*pull);
  }

#else
  for (int isys=0;isys < 2;isys++){
    pull = thePars->sysPar[isys]-1.;
    pull/=thePars->sysParUnc[isys];
    totL+=(0.5)*(pull*pull);
  }
  totL += thePars->cov->getLikelihood();
#endif

  ////////////////////////////////////////////////////////////////////////////////////////
  //Contribution from bias and smear priors
  // If using gaussion priors on bias and smear parameters, evelauate the likelihood here.
  // This is done by calling an atmfit pars method that will sum the contributions 
  if (flgUsePriorsInFit){
    totL += thePars->calcLogPriors();
  }



  return totL;
  
}
//*/

/*
////////////////////////////////////////////////
//Compute the total log liklihood by comparing all histograms
double histoCompare::getTotLnL(){

  double totL = 0.;


  ////////////////////////////////////////
  //contribution from histogram comparison  
  for (int isamp=0;isamp<nSamp;isamp++){
    for (int ibin=0;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
      	//    TH1D* hData = (TH1D*)hManager->getHistogramData(isamp,ibin,iatt)->Rebin(1,"hdata_rebinned");
      	//    TH1D* hPrediction = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt)->Rebin(1,"hmc_rebinned");
      	TH1D* hDataTmp = (TH1D*)hManager->getHistogramData(isamp,ibin,iatt);
      	//   TH1D* hPrediction = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt,0); //< get un-normalized histogram.
       	TH1D* hPrediction = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt,1); //< get normalized histogram.
	
	//    double hnorm = hDataTmp->Integral()/hPrediction->Integral();
	
	//    double partialL =  getLnL(hPrediction,hDataTmp,hnorm);
	//    cout<<"partialL :"<<partialL<<endl;     
	      double partialLnL = getLnL(hPrediction,hDataTmp);
	      totL+=partialLnL;
      }
    }
  }
  
  //////////////////////////////////////////////
  //contribution from flux/xsec priors
  double pull;
#ifndef T2K
  for (int isys=0;isys<thePars->nSysPars;isys++){
    // pull = thePars->sysPar[isys]-1.;
    // pull/=thePars->sysParUnc[isys];
    
    pull = thePars->getSysParameter(isys)-1.;
    pull/=thePars->sysParUnc[isys];
    
    totL+=(0.5)*(pull*pull);


  }

  // If using gaussion priors on bias and smear parameters, evelauate the likelihood here.
  // This is done by calling an atmfit pars method that will sum the contributions 
  if (flgUsePriorsInFit){
    totL += thePars->calcLogPriors();
  }

#else
  for (int isys=0;isys < 2;isys++){
    pull = thePars->sysPar[isys]-1.;
    pull/=thePars->sysParUnc[isys];
    totL+=(0.5)*(pull*pull);
  }
  totL += thePars->cov->getLikelihood();
#endif
  return totL;
  
}
*/

double histoCompare::getTotSumSq(){
  double totsumsq = 0.;
//  for (int isamp=0;isamp<nSamp;isamp++){
//    for (int ibin=0;ibin<nBin;ibin++){
//      for (int iatt=0;iatt<nAtt;iatt++){
//         //get modfied MC prediction
//         hMod = smearIt(hManager->hMC[isamp][ibin][0][iatt],Par[ibin][0][iatt][0],Par[ibin][0][iatt][1]);      
//         for (int icomp = 1;icomp<nComp;icomp++){
//           hMod->Add(smearIt(hManager->hMC[isamp][ibin][icomp][iatt],Par[ibin][icomp][iatt][0],Par[ibin][icomp][iatt][1]));
//         }
         //add error to total
        // hMod->Scale(Norm);
//         totsumsq+=getSumSq(hMod,hManager->hData[isamp][ibin][iatt]);
//      }
//    }
//  }
  return totsumsq;
}


double histoCompare::getSumSq(TH1D* h1, TH1D* h2){
  double sumsq = 0.;
  double diff;
  for (int ibin=10;ibin<=(h1->GetNbinsX()-10);ibin++){
    diff = h1->GetBinContent(ibin)-h2->GetBinContent(ibin);
    sumsq += (diff*diff);
  }
  return sumsq;
}


TH2D* histoCompare::show2DLnL(int parx, double xmin, double xmax, int pary, double ymin, double ymax, int npts){

  // what are the current par values? reset to these later
  double parxcurrent = thePars->getParameter(parx);
  double parycurrent = thePars->getParameter(pary);
  double LnLcurrent = getTotLnL();;

  // make a histogram
  TH2D* hL = new TH2D("hL","hL",npts,xmin,xmax,npts,ymin,ymax);

  // usefull for finding bin centers
  double xwidth = hL->GetXaxis()->GetBinWidth(1);
  double ywidth = hL->GetYaxis()->GetBinWidth(1);

  double minL = 1e6;
  double maxL = -1e5;

  // fill histo
  for (int ibin=1; ibin<npts; ibin++){
    for (int jbin=1; jbin<npts; jbin++){
      double xval = hL->GetXaxis()->GetBinLowEdge(ibin) + (xwidth/2.);
      double yval = hL->GetYaxis()->GetBinLowEdge(jbin) + (xwidth/2.);
      if (!thePars->fixPar[parx]){ 
       // cout<<"set par "<<parx<<"to "<<xval<<endl;
        thePars->setParameter(parx,xval);
      }
      if (!thePars->fixPar[pary]){
     //   cout<<"set par "<<pary<<"to "<<yval<<endl;
        thePars->setParameter(pary,yval);
      }
      double LnLval = getTotLnL();
   //   double diff = LnLval - LnLcurrent;
   //   if ((TMath::Abs(diff)>10) && (diff>0.)){ 
//        hL->SetBinContent(ibin,jbin,-1e5);
    //    continue;
     // }
      if (LnLval>maxL) maxL = LnLval;
      if (LnLval<minL) minL = LnLval;
      hL->SetBinContent(ibin,jbin,LnLval);
    }
  }
  cout<<"Min lnL: "<<minL<<endl;
  // fill diff 
  for (int ibin=1; ibin<npts; ibin++){
    for (int jbin=1; jbin<npts; jbin++){
      double binc = hL->GetBinContent(ibin,jbin);
      hL->SetBinContent(ibin,jbin,binc-minL);
    }
  }

  hL->GetZaxis()->SetRangeUser(0,5);
  hL->Draw("colz");
  
  thePars->setParameter(parx, parxcurrent);
  thePars->setParameter(pary, parycurrent);

  return hL;
}



TGraph2D* histoCompare::show2DLnLG(int parx, double xmin, double xmax, int pary, double ymin, double ymax, int npts){

  // what are the current par values? reset to these later
  float parxcurrent = thePars->getParameter(parx);
  float parycurrent = thePars->getParameter(pary);
  float LnLcurrent = getTotLnL();;

  // make a graph 
  const int NN = npts*npts;
  TGraph2D* g2 = new TGraph2D(NN);
  float X[NN];
  float Y[NN];
  float Z[NN];
  float dX = (xmax-xmin)/(float)npts;
  float dY = (ymax-ymin)/(float)npts;
  float xx = xmin;
  float yy = ymin;
  int nn = 0;
  float minL = 1e6;
  float maxL = -1e5;

  // fill arrays 
  for (int ibin=0; ibin<npts; ibin++){
    for (int jbin=0; jbin<npts; jbin++){
      Y[nn] = yy; 
      X[nn] = xx;
      if (!thePars->fixPar[parx]){ 
        thePars->setParameter(parx,X[nn]);
      }
      if (!thePars->fixPar[pary]){
        thePars->setParameter(pary,Y[nn]);
      }
      float LnLval = getTotLnL();
      if (LnLval>maxL) maxL = LnLval;
      if (LnLval<minL) minL = LnLval;
      Z[nn] = LnLval;
 //     g2->SetPoint(nn,xx,yy,LnLval);
      nn++;
      yy+=dY;
    }
    yy = ymin;
    xx+=dX;
  }
 
  // diff it
  for (int ip=0; ip<NN; ip++){
    g2->SetPoint(ip,X[ip],Y[ip],Z[ip] -minL);
  }

  // make it
//  TGraph2D* g2 = new TGraph2D(NN,X,Y,Z);

  cout<<"Min lnL: "<<minL<<endl;
  // fill diff 
  
  g2->GetZaxis()->SetRangeUser(0,5);
  g2->Draw("colz");
  
  thePars->setParameter(parx, parxcurrent);
  thePars->setParameter(pary, parycurrent);

  return g2;
}




////////////////////////////////////////////////////////////////////
//evalute log-likelihood between two histograms
double histoCompare::getLnL(TH1D* h1, TH1D* h2){
  double lnL = 0.;
//  double diff;
//  double term;
  double c1; //data
  double c2; //mc
  double errmc; //mcerr
  //double norm = hManager->normFactor; //normalization
  double norm = 1.0;
//  double dof=0.;
//  double quaderr;



  ///////////////////////////////////////////////////
  //assume poisson errors
 // for (int ibin=10;ibin<=(10);ibin++){
  int nedgebins= 5; //< need to ignore some edge bins
  for (int ibin=nedgebins;ibin<=(h1->GetNbinsX()-nedgebins);ibin++){
    c1 = h1->GetBinContent(ibin); //MC
    c2 = h2->GetBinContent(ibin); //data
    errmc = h1->GetBinError(ibin);

    if (c2<7) continue;
     
    lnL+=evalLnL(c2,c1,norm); //< tn186 likelihood definition
  //  lnL+=evalGausChi2WithError(c2,c1,errmc); //< tn186 likelihood definition

}
  return lnL;
}

void histoCompare::addHistogram(TH1D* h,int dataflg){
  if (dataflg==0){
    hMC[nMCHist] = h;
    nMCHist++;
  }
  else {
    hData[nDataHist] = h;
    nDataHist++;  
  }  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Read in histograms from a specified file
void histoCompare::readFromFile(const char* filerootname,int nsamp, int nbin, int ncomp, int natt){

  // set variables
  nSamp = nsamp;
  nBin  = nbin;
  nComp = ncomp; 
  nAtt  = natt;

  // read in histograms
  hManager = new histoManager(filerootname,nsamp,nbin,ncomp,natt);
  

  /*
  double ndataevents=0;
  double nmcevents=0;
  double events;
  //count total events
  for (int jsamp=0;jsamp<nSamp;jsamp++){
    for (int jbin=0;jbin<nBin;jbin++){
      for (int jatt=0;jatt<nAtt;jatt++){
        events =  hManager->hData[jsamp][jbin][jatt]->Integral();
        cout<<"histo "<<jsamp<<"-"<<nbin<<"-"<<jatt<<" has "<<events<<" events."<<endl;
        ndataevents+=events;
        for (int jcomp=0;jcomp<nComp;jcomp++){
          events = hManager->hMC[jsamp][jbin][jcomp][jatt]->Integral();
          cout<<"MC histo "<<jsamp<<"-"<<jbin<<"-"<<jcomp<<"-"<<jatt<<" has "<<events<<" events."<<endl;
          nmcevents+=events;
        }
      }
    }
  }


  Norm = ndataevents/nmcevents;;
  */

  return;
}


void histoCompare::readFromFile(const char* filerootname,int nsamp, int nbin, int ncomp, int nmode, int natt){
  nSamp = nsamp;
  nBin  = nbin;
  nComp = ncomp; 
  nAtt  = natt;
  nMode = nmode;
  hManager = new histoManager(filerootname,nsamp,nbin,ncomp,natt,nmode,true);
  double ndataevents=0;
  double nmcevents=0;
  double events;
}


void histoCompare::makeFakeData(){
   // sets all modified histograms as data 
   for (int isamp=0; isamp<nSamp; isamp++){
     for (int ibin=0; ibin<nBin; ibin++){
       for (int iatt=0; iatt<nAtt; iatt++){
         TString hname = hManager->hData[isamp][ibin][iatt]->GetName();
         hManager->hData[isamp][ibin][iatt] = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt,0)
         ->Clone(hname.Data());
       }
     }
   }
   hManager->normFactor = 1.0;
   return;
}

/////////////////////////////////////////////////////////////////////
//initializes all necessary compoents
void  histoCompare::initialize(histoManager* hm, atmFitPars* apars){
  cout<<"histoCompare: Initialization: "<<endl;
  thePars = apars;
  hManager = hm;
  nComp = thePars->nComponents;
  nBin  = thePars->nBins;
  nAtt  = thePars->nAttributes;
  nSamp = thePars->nSamples;
  cout<<"    bins: "<<nBin<<endl;
  cout<<"    sampless: "<<nSamp<<endl;
  cout<<"    components: "<<nComp<<endl;
  cout<<"    attributes: "<<nAtt<<endl; 
  return; 
}

void histoCompare::setupPars(int nsyspars){
  thePars = new atmFitPars(nSamp,nBin,nComp,nAtt,"tn186");
 // thePars->setNorm(Norm);
  hManager->setFitPars(thePars);
  return;
}

void histoCompare::setupPars(atmFitPars *a)
{
  thePars = a;
  hManager->setFitPars(thePars);
}

histoCompare::histoCompare(){
  nameTag = "manualSetup";
  nMCHist=0;
  nDataHist=0;
  cout<<"created comparison object: "<<nameTag.Data()<<endl;
  //setup initial debug
  int jhist = 0;
  cc = new TCanvas("cc","cc",700,600);
  while (jhist<10){
      parDebug[jhist][0] = 1.0;
      parDebug[jhist][1] = 0.0;
      jhist++;
  }
  useLnLType=0;
  return;
}


//////////////////////////////////////////////////////////
// Use this constructor 
histoCompare::histoCompare(const char* parfile, bool sep)
  : separateNeutMode(sep)
{

  //read in parameter file
  runPars = new sharedPars(parfile);
  runPars->readParsFromFile();
  nameTag= runPars->globalRootName.Data();

  //setup canvas
  cc = new TCanvas("cc","cc",700,600);
  
  //set LnL definition to default (0)
  useLnLType=0;

  //MCMC tuning parameter
  tunePar = runPars->MCMCTunePar;

  //Use smear parameters or no?
  flgFixAllSmearPars = runPars->flgFixAllSmearPars;

  //Should we use priors in this fit?
  flgUsePriorsInFit = runPars->flgUsePriorsInFit;

  //MCMC nsteps;
  MCMCNSteps = runPars->MCMCNSteps;

  //read in pre-filled histograms using histoManager
  int nbins = runPars->nFVBins;
  int ncomponents = runPars->nComponents;
  int nsamples = runPars->nSamples;
  int nattributes  = runPars->nAttributes;
  TString histofilename = runPars->hFactoryOutput;
  readFromFile(histofilename.Data(),nsamples,nbins,ncomponents,nattributes); 
  
  //setup fit parameters
  thePars = new atmFitPars(parfile);
  hManager->setFitPars(thePars);
  int nsyspars = thePars->nSysPars;
  if (flgFixAllSmearPars){
    thePars->fixAllSmearPars();
  }

  
  //read in splines if you're into that
  if (runPars->useSplinesFlg){
    setupSplines(runPars->splineFactoryOutput.Data());
  };


}


void histoCompare::useFakeData(){
 return; 
}















#endif
