#ifndef HISTOCOMPARE_C
#define HISTOCOMPARE_C

#include "histoCompare.h"
#include "markovTools.C"
#include <time.h>

histoCompare* histoCompare::staticthis;

void histoCompare::readFitPars(const char* filename){
  thePars->readPars(filename);
  return;
}

void histoCompare::timetest(int ntry){
  clock_t t1,t2;
  float par[100];
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
  float diff;
  float result;
  int itry=0;
  t1=clock();
  while (itry<ntry){
    //getTotLnL1D(result,par); 
    getTotLnL();
    itry++;
  }
  t2=clock();
  diff = ((float)t2-(float)t1)/(float)ntry;
  float diff1 = diff;
  cout<<"time 1: "<<diff<<endl;
  itry=0;
  t1=clock();
  while (itry<ntry){
    //getTotLnL1D(result,par); 
    getTotLnL();
    itry++;
  }
  t2=clock();
  diff = ((float)t2-(float)t1)/(float)ntry;
  cout<<"time 2: "<<diff<<endl;
  cout<<"changes: "<<diff1-diff<<endl;
  return;
}

void histoCompare::tuneMCMC(int ncycles,int nsteps,float goal){

  //setup mc mcmctools
//  tunePar = guess;
  const int npars = thePars->nTotPars; //< total number of parameters in fit
  float par[npars]; //< container for parameters
  int parindex = 0;
  float result = 0.;
  markovTools* mc = new markovTools(npars); //< create markovTools object
  mc->setTuneParameter(tunePar);
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<npars;ipar++){
    par[ipar]=thePars->getParameter(ipar);
    mc->setParVar(ipar,thePars->parUnc[ipar]);
  }
  //set initial state
  getTotLnL1D(result,npars, par);//< get total likelihood from 1D array
  mc->setL((float)result);//< sets the initial likelihood
  float Linit = result;
  //run tuning
  for (int icycle=0;icycle<ncycles;icycle++){
    float xaccepted=0.;
    int   istep=0;
    while (istep<nsteps){
      istep = mc->iStep;
      mc->proposeStep(par); //< propose a new step
      getTotLnL1D(result, npars,par);  //< get likelihood of new step
      cout<<result-Linit<<endl;
      if (mc->acceptStepLnL(result,par)){ //< check if new step is accepted
        xaccepted++; 
      }
    }
    float rate = xaccepted/(float)nsteps;
    cout<<"acceptance rate: "<<rate<<endl;
    tunePar*=(rate/goal);
    cout<<"new tune parameter: "<<tunePar<<endl;
  }
  return; 
}

void histoCompare::runMCMC(int nsteps){
  //setup mcmc tools
  const int npars = thePars->nTotPars; //< total number of parameters in fit
  float par[npars]; //< container for parameters
  int parindex = 0;
  float result = 0.;
  float tune = 0.11;
  markovTools* mc = new markovTools(npars); //< create markovTools object
  mc->setTuneParameter(tune);
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<npars;ipar++){
    par[ipar]=thePars->getParameter(ipar);
    mc->setParVar(ipar,(errParHi[ipar]-errParLo[ipar]));
  }
  //set initial state
  getTotLnL1D(result,npars, par);//< get total likelihood from 1D array
  mc->setL((float)result);//< sets the initial likelihood
  int currentstep=0;
  while (currentstep<nsteps){
    currentstep = mc->iStep;
    mc->proposeStep(par); //< fills par[] with proposed params
    getTotLnL1D(result, npars,par);   
    mc->acceptStepLnL(result,par); //< if step is accepted, istep++, and written to histo
  }
  mc->savePath("mcmctest.root");
  return;
}

/*
float histoCompare::getErrHi(int ibin,int icomp,int iatt,int imod){
  if (fixPar[ibin][icomp][iatt][imod]==1) return 0.;
  float thresh = 1.0;
  float Ldiff = 0.;
  float Lbest = getTotLnL();
  float parbest = Par[ibin][icomp][iatt][imod];
  float parval = 0.;
  float dpar = 1.;
  float loerr;
  int ntry = 0;
  int ntrymax=1000;
  if (imod==0) dpar = 0.005;
  if (imod==1) dpar = 10.;
  //course search
  while ((Ldiff<1)&&(ntry<ntrymax)){
    Par[ibin][icomp][iatt][imod]+=dpar; //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
  //  cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
  }
 // cout<<"ntry: "<<ntry<<endl;
  Par[ibin][icomp][iatt][imod]-=dpar;
  Ldiff = 0;
  dpar*=0.2;
  ntry=0;
  while ((Ldiff<1)&&(ntry<ntrymax)){
    Par[ibin][icomp][iatt][imod]+=dpar; //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
  //  cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
  }
  Par[ibin][icomp][iatt][imod]-=dpar;
  Ldiff = 0;
  dpar*=0.1;
  ntry=0;
  while ((Ldiff<1)&&(ntry<ntrymax)){
    Par[ibin][icomp][iatt][imod]+=dpar; //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
  //  cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
  }

  loerr = Par[ibin][icomp][iatt][imod]-parbest;
  Par[ibin][icomp][iatt][imod] = parbest;
  return loerr; 
}
*/

void histoCompare::saveFitPars(const char* filename){
  thePars->savePars(filename);
  return;
}

float histoCompare::getErrHi(int ipar){
  //scan through log likelihood by decreasing parameter until threshold is reached
  float thresh = 1.0; //likelihood threshold
  float Ldiff = 0.;  //difference in likelihood from current value
  float Lbest = getTotLnL(); //current likelihood value
  float parbest = thePars->pars[ipar];
  float parval = thePars->pars[ipar]; 
  float dpar = thePars->parUnc[ipar]/2.;
  float hierr;
  int ntry = 0;
  int ntrymax=1000;
  //coarse search
  while (Ldiff<1){
    parval+=dpar;
    thePars->setParameter(ipar,parval); //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
    ntry++;
  }
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
 // cout<<"ntry: "<<ntry<<endl;
  hierr = thePars->pars[ipar];
  thePars->setParameter(ipar,parbest);
  return hierr; 
}



float histoCompare::getErrLo(int ipar){
  //scan through log likelihood by decreasing parameter until threshold is reached
  float thresh = 1.0; //likelihood threshold
  float Ldiff = 0.;  //difference in likelihood from current value
  float Lbest = getTotLnL(); //current likelihood value
  float parbest = thePars->pars[ipar];
  float parval = thePars->pars[ipar]; 
  float dpar = thePars->parUnc[ipar]/10.;
//  cout<<"dpar: "<<dpar<<endl;
  float loerr;
  int ntry = 0;
  int ntrymax=1000;
  //coarse search
  while (Ldiff<1){
    parval-=dpar;
    thePars->setParameter(ipar,parval); //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
 //   cout<<"ntry: "<<ntry<<endl;
 //   cout<<"Ldiff: "<<Ldiff<<endl;
 //   cout<<"par: "<<thePars->pars[ipar]<<endl;
    ntry++;
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
  }
//  cout<<"ntry: "<<ntry<<endl;
  parval+=dpar;
  thePars->setParameter(ipar,parval);
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  //very fine search
  while ((Ldiff<1)&&(ntry<ntrymax)){
    parval-=dpar;
    thePars->setParameter(ipar,parval); //modify paramete
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
    ntry++;
  }
//  cout<<"ntry: "<<ntry<<endl;
  loerr = thePars->pars[ipar];
  thePars->setParameter(ipar,parbest);
  return loerr; 
}

/*
float histoCompare::getErrLo(int ibin,int icomp,int iatt,int imod){
  //scan through log likelihood by decreasing parameter
  if (fixPar[ibin][icomp][iatt][imod]==1) return 0.;
  float thresh = 1.0;
  float Ldiff = 0.;
  float Lbest = getTotLnL();
  float parbest = Par[ibin][icomp][iatt][imod];
  float parval = 0.;
  float dpar = 1.;
  float loerr;
  int ntry = 0;
  int ntrymax=1000;
  if (imod==0) dpar = 0.005;
  if (imod==1) dpar = 10.;
  //course search
  while (Ldiff<1){
    Par[ibin][icomp][iatt][imod]-=dpar; //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
  //  cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
  }
 // cout<<"ntry: "<<ntry<<endl;
  Par[ibin][icomp][iatt][imod]+=dpar;
  Ldiff = 0;
  dpar*=0.2;
  ntry=0;
  while ((Ldiff<1)&&(ntry<ntrymax)){
    Par[ibin][icomp][iatt][imod]-=dpar; //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
  //  cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
  }
  Par[ibin][icomp][iatt][imod]+=dpar;
  Ldiff = 0;
  dpar*=0.1;
  ntry=0;
  while ((Ldiff<1)&&(ntry<ntrymax)){
    Par[ibin][icomp][iatt][imod]-=dpar; //modify parameter
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
  //  cout<<"Ldiff: "<<Ldiff<<endl;
    ntry++;
  }
  loerr = parbest-Par[ibin][icomp][iatt][imod];
  Par[ibin][icomp][iatt][imod] = parbest;
  return loerr; 
}
*/

void histoCompare::profileL(int ipar, float range, int npts){
  TString pname = "p";
  pname.Append("_profile.png");
  float bestpoint = thePars->pars[ipar];
  cout<<"Fitted value: "<<bestpoint<<endl;
  float dx = range/(float)npts;
  float xx = bestpoint - (range/2.);
  float ll;
  float lbest = getTotLnL();
  if (hProf) hProf->Delete();
  hProf = new TH1F("hprof","hprof",npts,xx,(xx+range));
  for (int ipoint=0;ipoint<npts;ipoint++){
   // cout<<"filling point" <<ipoint<<endl;
    //Par[ibin][icomp][iatt][imod] = xx;
    thePars->setParameter(ipar,xx);
    ll = getTotLnL();
    hProf->SetBinContent(ipoint+1,ll-lbest);
//    hProf->SetBinContent(ipoint+1,ll);

    xx+=dx;
  } 
  //set parameter back to initial value
  thePars->setParameter(ipar,bestpoint);
  hProf->SetLineColor(9);
  hProf->GetYaxis()->SetTitle("#Delta #chi^{2}");
  TString xname = "parameter ";
  xname.Append(Form("%d",ipar));
  hProf->GetXaxis()->SetTitle(xname.Data());
  hProf->Draw("c");
  cc->Print(pname.Data());
  return;
}




void histoCompare::profileL(int ibin, int icomp, int iatt, int imod, float range, int npts){
  TString pname = "p";
  pname.Append("_profile.png");
  //float bestpoint = Par[ibin][icomp][iatt][imod];
  float bestpoint = thePars->histoPar[ibin][icomp][iatt][imod];

  cout<<"Fitted value: "<<bestpoint<<endl;
  float dx = range/(float)npts;
  float xx = bestpoint - (range/2.);
  float ll;
  float lbest = getTotLnL();
  if (hProf) hProf->Delete();
  hProf = new TH1F("hprof","hprof",npts,xx,(xx+range));
  for (int ipoint=0;ipoint<npts;ipoint++){
   // cout<<"filling point" <<ipoint<<endl;
    //Par[ibin][icomp][iatt][imod] = xx;
    thePars->setParameter(ibin,icomp,iatt,imod,xx);
    ll = getTotLnL();
    hProf->SetBinContent(ipoint+1,ll-lbest);
 //   hProf->SetBinContent(ipoint+1,ll);

    xx+=dx;
  } 
  //set parameter back to initial value
  thePars->setParameter(ibin,icomp,iatt,imod,bestpoint);
  hProf->SetLineColor(9);
  hProf->Draw("c");
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
  hPar  = new TH1F("hpar","hpar",nbinstot,0,nbinstot);
  hParErrLo = new TH1F("hparerrlo","hparerrlo",nbinstot,0,nbinstot);
  hParErrHi = new TH1F("hparerrhi","hparerrhi",nbinstot,0,nbinstot);
  for (int icomp=0;icomp<nComp;icomp++){
    parindex=thePars->getParIndex(ibin,icomp,iatt,imod);
    hPar->SetBinContent(icomp+1,thePars->getParameter(parindex));
    X[icomp]=(double)hPar->GetBinCenter(icomp+1);
    Y[icomp]=(double)thePars->getParameter(parindex);
 //   hPar->GetXaxis()->SetBinLabel(icomp+1,parName[ibin][icomp][iatt][imod].Data());
    EYL[icomp]=(double)thePars->getParameter(parindex) -errParLo[parindex];
    EYH[icomp]=errParHi[parindex] - (double)thePars->getParameter(parindex);
    hParErrLo->SetBinContent(icomp+1,(float)EYL[icomp]);
    hParErrHi->SetBinContent(icomp+1,(float)EYH[icomp]);
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

void histoCompare::showModHiso(int isamp,int ibin, int icomp, int iatt, float smear, float bias){
  float biastmp = thePars->histoPar[ibin][icomp][iatt][1];
  float smeartmp = thePars->histoPar[ibin][icomp][iatt][1];
  thePars->setParameter(ibin,icomp,iatt,0,smear);
  thePars->setParameter(ibin,icomp,iatt,1,bias);
  hMod = hManager->getModHistogram(isamp,ibin,icomp,iatt); //gets the modified histogram
  hTmp = hManager->getHistogram(isamp,ibin,icomp,iatt);
//  if (hMod) hMod->Delete();
//  if (hTmp) hTmp->Delete();
//  if (hTot) hTot->Delete();
//  hMod = (TH1F*) hManager->hMC[isamp][ibin][icomp][iatt]->Clone("hclonetmp");
//  hTot = (TH1F*) hManager->hMC[isamp][ibin][icomp][iatt]->Clone("htottmp");
//  hTmp = (TH1F*) hManager->hMC[isamp][ibin][icomp][iatt]->Clone("htmp");
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
  hMod->Draw("same");
  //set parameter back to original
  thePars->setParameter(ibin,icomp,iatt,0,smeartmp);
  thePars->setParameter(ibin,icomp,iatt,1,biastmp);

  return;
}

TH1F* histoCompare::showSmear(TH1F* h, float smear, float bias){
  TH1F* hh = smearIt(h,smear,bias);
  return hh;
}

void histoCompare::showFitResult(int isamp,int ibin,int iatt){
  //sum up mc components
  //get MC prediction
  hMod = (TH1F*)hManager->getSumHistogramMod(isamp,ibin,iatt)->Clone("hmod");
  hTmp = hManager->getSumHistogram(isamp,ibin,iatt);
 // hMod = smearIt(hManager->hMC[isamp][ibin][0][iatt],smear,bias);
//  hTot = (TH1F*)hManager->hMC[isamp][ibin][0][iatt]->Clone("htot");
//  for (int jcomp=1;jcomp<nComp;jcomp++){
//    hTot->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
//    smear = bestPar[ibin][jcomp][iatt][0];
//    bias  = bestPar[ibin][jcomp][iatt][1];
 //   hMod->Add(smearIt(hManager->hMC[isamp][ibin][jcomp][iatt],smear,bias));
//  }
//  hTot->SetLineColor(kRed);
 // hMod->Scale(Norm);
 // hTmp->Scale(Norm);
  hMod->SetLineColor(kRed);
  hTmp->SetLineColor(kBlue);
//  hTot->Rebin(rebinFactor);
//  hMod->Rebin(rebinFactor);
  hMod->Draw();
  hTmp->Draw("same");

  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->Draw("samee");
  return;
}

void histoCompare::showFitEffect(int isamp,int ibin,int icomp,int iatt){
  //sum up mc components
  hMod = hManager->getModHistogram(isamp,ibin,icomp,iatt); //gets the modified histogram
  hTmp = hManager->getSumHistogram(isamp,ibin,iatt);

//  float smear = bestPar[ibin][icomp][iatt][0];
//  float bias  = bestPar[ibin][icomp][iatt][1];
//  cout<<"SMEAR: "<<smear<<endl;
//  cout<<"BIAS:  "<<bias<<endl;
//  hMod = smearIt(hManager->hMC[isamp][ibin][icomp][iatt],smear,bias);
//  hTot = (TH1F*)hManager->hMC[isamp][ibin][icomp][iatt]->Clone("htot");
  for (int jcomp=0;jcomp<nComp;jcomp++){
    if (jcomp!=icomp){
      hMod->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
    }
  }
  hMod->SetLineColor(kRed);
 // hTmp->Scale(Norm);
  hTmp->SetLineColor(kBlue);
//  hMod->Scale(Norm);
//  hTot->Rebin(rebinFactor);
//  hMod->Rebin(rebinFactor);
//  hTmp->Draw();
  hMod->Draw();
  hTmp->Draw("same");
  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->Draw("samee");
  return;
}

void histoCompare::showFitHisto(int isamp,int ibin,int icomp,int iatt){
  float smear = bestPar[ibin][icomp][iatt][0];
  float bias  = bestPar[ibin][icomp][iatt][1];
  cout<<"SMEAR: "<<smear<<endl;
  cout<<"BIAS:  "<<bias<<endl;
  hMod = smearIt(hManager->hMC[isamp][ibin][icomp][iatt],smear,bias);
  hManager->hMC[isamp][ibin][icomp][iatt]->SetLineColor(kRed);
  hManager->hMC[isamp][ibin][icomp][iatt]->Draw();
  hMod->SetLineColor(kBlue);
//  hMod->Smooth(1);
//  convolveThisHisto(*hMod,hMod->GetBinWidth(2)*0.5,0.);
  hMod->Draw("same");
  return;
}

void histoCompare::LnLPreFit(){

  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //threshold to determine if a parameter is fit or not
  //if the MC histograms for this parameter have a size less than this value,
  //don't bother fitting them!
  float nthresh = 100.;

  //sets the precision of the fits
  double parerr = 0.005;  
  
  //individually fit each parameter
  int parindex =0;
  int parindextmp;
  //parameter name container
  TString parnametmp;  

  float  parinit; //container to temporarily store initial values 

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
        hTot = (TH1F*)hManager->hMC[0][ibin][icomp][iatt]->Clone("htot");
        for (int isamp=1;isamp<nSamp;isamp++){
          hTot->Add(hManager->hMC[isamp][ibin][icomp][iatt]);
        }
        cout<<"total entries: "<<hTot->GetEntries()<<endl;
        //check to make sure histogram is above fitting threshold
        if (hTot->GetEntries()<nthresh){
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][0].Data()<<" (ENTRIES TOO LOW!) "<<endl; 
          fixPar[ibin][icomp][iatt][0]=1;
          thePars->fixParameter(ibin,icomp,iatt,0);
          thePars->fixParameter(ibin,icomp,iatt,1);
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][1].Data()<<" (ENTRIES TOO LOW!) "<<endl; 
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

//  parindex = 0;
//  for (int kbin=0;kbin<nBin;kbin++){
//    for (int katt=0;katt<nAtt;katt++){
//      for (int kcomp=0;kcomp<nComp;kcomp++){
//        fit->SetParameter(parindex,parName[kbin][kcomp][katt][0].Data(),Par[kbin][kcomp][katt][0],parerr,0,0);
//        parindex++;
//        bestPar[kbin][kcomp][katt][1] = Par[kbin][kcomp][katt][1];
//        fit->SetParameter(parindex,parName[kbin][kcomp][katt][1].Data(),Par[kbin][kcomp][katt][1],parerr,0,0);
//        parindex++;
//      }
//    }
//  }

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }

//  debugging
//  fit->ReleaseParameter(npars-1);
//  fit->ExecuteCommand("MIGRAD",0,0);
//  fit->FixParameter(npars-1);
//  fit->ExecuteCommand("MIGRAD",0,0);


  //fit flux ans xsec
  parindex = thePars->nTotPars-thePars->nSysPars;
  for (int isyspar=0;isyspar<thePars->nSysPars;isyspar++){
    fit->ReleaseParameter(parindex);
    parindex++;
  }
  fit->ExecuteCommand("SIMPLEX",0,0);
  parindex = 0;
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
 
  //run individual bias fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
          parindex++; //starts on odd parameter (bias only)
          //release biasparameter to be fit
          if (thePars->checkFixFlg(jbin,jcomp,jatt,1)==1){
            parindex++; //do nothing if parameter is fixed
            continue;
          }
          fit->ReleaseParameter(parindex);
          cout<<"fitting "<<jbin<<jcomp<<jatt<<1<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,1)<<endl;
          fit->ExecuteCommand("SIMPLEX",0,0);
          fit->FixParameter(parindex);
       //   bestPar[jbin][jcomp][jatt][1] = fit->GetParameter(parindex); 
          parindex++;
      }
    }
  }
  parindex=0;
  //run individual smear fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
          //starts on even parameter (bias only)
          //release smear parameter to be fit
          if (thePars->checkFixFlg(jbin,jcomp,jatt,0)==1){
            parindex++; //do nothing if parameter is fixed
            continue;
          }
          fit->ReleaseParameter(parindex);
          cout<<"fitting "<<jbin<<jcomp<<jatt<<0<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,0)<<endl;
          fit->ExecuteCommand("SIMPLEX",0,0);
          fit->FixParameter(parindex);
         // bestPar[jbin][jcomp][jatt][0] = fit->GetParameter(parindex); 
          parindex+=2;
      }
    }
  }


  //print final results
  for (int ipar=0;ipar<npars;ipar++){
    thePars->setParameter(ipar,fit->GetParameter(ipar));
    cout<<"  PAR "<<ipar<<" FIT RESULT: "<<thePars->pars[ipar]<<endl;
  }



/*
  for (int pbin=0;pbin<nBin;pbin++){
    for (int patt=0;patt<nAtt;patt++){
      for (int pcomp=0;pcomp<nComp;pcomp++){
        for (int pmod=0;pmod<2;pmod++){
          cout<<"  PAR "<<parName[pbin][pcomp][patt][pmod].Data()<<" FIT RESULT: ";
          cout<<thePars->histoPar[pbin][pcomp][patt][pmod]<<endl;
        }
      }
    }
  }
  for (int isys=0;isys<thePars->nSysPars;isys++){
          cout<<thePars->sysPar[isys]<<endl;
  }
 */

  cout<<"$$$$$$$$$$$$$ END LNL PRE FIT $$$$$$$$$$$$$"<<endl;

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
  fit->ExecuteCommand("MIGRAD",0,0); 
  //print results
  for (int jpar=(thePars->nTotPars-thePars->nSysPars);jpar<thePars->nTotPars;jpar++){
    cout<<"PAR: "<<jpar<<" "<<fit->GetParameter(jpar)<<endl;
  }

  return;
}



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
  fit->ReleaseParameter(ipar);
  //fit that thang
  fit->ExecuteCommand("SIMPLEX",0,0); 
  fit->ExecuteCommand("MIGRAD",0,0); 
  //print results
  cout<<"RESULT: "<<fit->GetParameter(ipar)<<endl;
 
  return;
}

void histoCompare::LnLFit(){
  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //threshold to determine if a parameter is fit or not
  //if the MC histograms for this parameter have a size less than this value,
  //don't bother fitting them!
  float nthresh = 50.;

  //sets the precision of the fits
  double parerr = 0.001;  
  
  //individually fit each parameter
  int parindex =0;

  //parameter name container
  TString parnametmp;  

  float  parinit; //container to temporarily store initial values 

  //total number of parameters to be fit
  int npars = thePars->nTotPars;

  //run the prefit
  LnLPreFit();

  cout<<"$$$$$$$$$$$$$$$$$ LNL FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //setup the fitter!
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
    int kbin=thePars->binOfPar[ipar];
    int kcomp=thePars->compOfPar[ipar];
    int katt = thePars->attOfPar[ipar];
    aname = "parameter_";
    aname.Append(ipar);
    fit->SetParameter(ipar,aname.Data(),thePars->pars[ipar],parerr,0,0);
  }
 
  parindex = 0;
  int parindextmp = 0;
  //do individual fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      //start of fit block//

      //fix all parameters
      for (int jpar=0;jpar<npars;jpar++){
        fit->FixParameter(jpar);
      }

      //for bin 0, float the xsec and flux parameters as well//
      if (jbin==0){
        //release all flux and xsec pars
        for (int isyspar=(thePars->nTotPars-thePars->nSysPars);isyspar<thePars->nSysPars;isyspar++){
          fit->ReleaseParameter(isyspar);
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



/*
  //print final results
  cout<<"  ----------------------------------------  "<<endl;
  for (int pbin=0;pbin<nBin;pbin++){
    for (int patt=0;patt<nAtt;patt++){
      for (int pcomp=0;pcomp<nComp;pcomp++){
        for (int pmod=0;pmod<2;pmod++){
          cout<<"  PAR "<<parName[pbin][pcomp][patt][pmod].Data()<<" FIT RESULT: ";
          cout<<thePars->histoPar[pbin][pcomp][patt][pmod]<<endl;
     //     errParLo[pbin][pcomp][patt][pmod] = getErrLo(pbin,pcomp,patt,pmod);
     //     errParHi[pbin][pcomp][patt][pmod] = getErrHi(pbin,pcomp,patt,pmod);
        }
      }
    }
  }
*/

  //print final results
  for (int ipar=0;ipar<npars;ipar++){
    thePars->setParameter(ipar,fit->GetParameter(ipar));
    errParLo[ipar]=getErrLo(ipar);
    errParHi[ipar]=getErrHi(ipar);
    thePars->parUnc[ipar]=(errParHi[ipar]-errParLo[ipar]);
    cout<<"  PAR "<<ipar<<" FIT RESULT: "<<thePars->pars[ipar]<<" +/- : "<<thePars->parUnc[ipar]<<endl;
  }

  //end
  cout<<"  ----------------------------------------  "<<endl;
  cout<<"$$$$$$$$        FIT  COMPLLETE      $$$$$$$$"<<endl;

  return;


}



void histoCompare::sumSqPrefit(){
  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //threshold to determine if a parameter is fit or not
  //if the MC histograms for this parameter have a size less than this value,
  //don't bother fitting them!
  float nthresh = 50.;

  //sets the precision of the fits
  double parerr = 0.05;  
  
  //individually fit each parameter
  int parindex =0;

  //parameter name container
  TString parnametmp;  

  float  parinit; //container to temporarily store initial values 

  //total number of parameters to be fit
  int npars = nBin*nComp*nAtt*2;

  cout<<"$$$$$$$$ CALLED SUM SQUARES PREFIT $$$$$$$$"<<endl;
  cout<<"  ----------------------------------------  "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ----------------------------------------  "<<endl;


  //fix parameters with too few events to be fit
  parindex = 0;
  for (int ibin=0;ibin<nBin;ibin++){
    for (int icomp=0;icomp<nComp;icomp++){
      for (int iatt=0;iatt<nAtt;iatt++){
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
        hTot = (TH1F*)hManager->hMC[0][ibin][icomp][iatt]->Clone("htot");
        for (int isamp=0;isamp<nSamp;isamp++){
          hTot->Add(hManager->hMC[isamp][ibin][icomp][iatt]);
        }
        //check to make sure histogram is above fitting threshold
        if (hTot->GetEntries()<nthresh){
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][0].Data()<<" (ENTRIES TOO LOW!) "<<endl; 
          fixPar[ibin][icomp][iatt][0]=1;
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][1].Data()<<" (ENTRIES TOO LOW!) "<<endl; 
          fixPar[ibin][icomp][iatt][1]=1;
        }
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
  parindex = 0;
  for (int kbin=0;kbin<nBin;kbin++){
    for (int kcomp=0;kcomp<nComp;kcomp++){
      for (int katt=0;katt<nAtt;katt++){
        fit->SetParameter(parindex,parName[kbin][kcomp][katt][0].Data(),Par[kbin][kcomp][katt][0],parerr,0,0);
        parindex++;
        fit->SetParameter(parindex,parName[kbin][kcomp][katt][1].Data(),Par[kbin][kcomp][katt][1],parerr,0,0);
        parindex++;
      }
    }
  }
  
  parindex = 0;
  //do individual fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jcomp=0;jcomp<nComp;jcomp++){
      for (int jatt=0;jatt<nAtt;jatt++){
        for (int jmod=0;jmod<2;jmod++){

          //fill initial value to restore parameter array later
          parinit = Par[jbin][jcomp][jatt][jmod];

          //if parameter was previously fixed, don't bother fitting it
          if (fixPar[jbin][jcomp][jatt][jmod]==1){
            bestPar[jbin][jcomp][jatt][jmod] = parinit;
            parindex++;
            continue;
          }

          //fit a single parameter
          fit->ReleaseParameter(parindex);   
          //fix all other parameters
          for (int jpar=0;jpar<npars;jpar++){
             if (jpar!=parindex) fit->FixParameter(jpar); 
          }
          fit->ExecuteCommand("SIMPLEX",0,0); //run the fit
   //       fit->ExecuteCommand("MIGRAD",0,0); //run the fit

          //get fit result
          bestPar[jbin][jcomp][jatt][jmod]=fit->GetParameter(parindex);

          //restore default parameters
          Par[jbin][jcomp][jatt][jmod] = parinit;
          fit->SetParameter(parindex,parName[jbin][jcomp][jatt][jmod].Data(),Par[jbin][jcomp][jatt][jmod],parerr,0,0);



          //running count of 1D parameter index
          parindex++;
        }
      }
    }
  }

 
  for (int pbin=0;pbin<nBin;pbin++){
    for (int pcomp=0;pcomp<nComp;pcomp++){
      for (int patt=0;patt<nAtt;patt++){
        for (int pmod=0;pmod<2;pmod++){
          cout<<"  PAR "<<parName[pbin][pcomp][patt][pmod].Data()<<" FIT RESULT: ";
          cout<<Par[pbin][pcomp][patt][pmod]<<" -> "<<bestPar[pbin][pcomp][patt][pmod]<<endl;
        }
      }
    }
  }
  
  cout<<"  ----------------------------------------  "<<endl;
  cout<<"$$$$$$$$        FIT  COMPLLETE      $$$$$$$$"<<endl;

  return;
}


float histoCompare::getNDiff(){
  hTot = smearIt(hMC[0],parDebug[0][0],parDebug[0][1]);
  for (int i=1;i<nMCHist;i++){
    hTot->Add(smearIt(hMC[i],parDebug[i][0],parDebug[i][1]));
  }
  float xcut = 100.;
  int   cutbin = hTot->FindBin(xcut);
  int   binmax = hTot->GetNbinsX();
  float diff1 = hTot->Integral(1,cutbin)-hData[0]->Integral(1,cutbin);
  float diff2 = hTot->Integral(cutbin,binmax)-hData[0]->Integral(cutbin,binmax);
  return diff1*diff1 + diff2*diff2;
}

void histoCompare::showMod(int imchist){
  hMod = smearIt(hMC[imchist],parDebug[imchist][0],parDebug[imchist][1]);
  hMod->Draw();
  return;
}

void histoCompare::drawResult(int ihist){
  hModDebug = smearIt(hMC[ihist],parDebug[ihist][0],parDebug[ihist][1]);
  hModDebug->SetLineColor(kBlue);
  hModDebug->Draw();
  hMC[ihist]->SetLineColor(kRed);
  hMC[ihist]->Draw("same");
  hData[ihist]->Draw("same");
  return;
}

void histoCompare::getTotLnL1D(float& result,int npar, float par[]){

  for (int ipar=0;ipar<npar;ipar++){
    thePars->setParameter(ipar,par[ipar]);
  }

//  int index=0;
//  for (int ibin=0;ibin<nBin;ibin++){
 //   for (int icomp=0;icomp<nComp;icomp++){
 //     for (int iatt=0;iatt<nAtt;iatt++){
  //      Par[ibin][icomp][iatt][0] = par[index];
  //      Par[ibin][icomp][iatt][1] = par[index+1]; 
  //      index+=2;
  //    }
  //  }
 // }
 // */
  result = getTotLnL();
}


void histoCompare::lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg){

  for (int ipar=0;ipar<staticthis->thePars->nTotPars;ipar++){
    staticthis->thePars->setParameter(ipar,par[ipar]);
  }

 // set MD parameter array to match the values of the internal pars[] of the fittter
//  int index=0;
//  for (int ibin=0;ibin<staticthis->nBin;ibin++){
//    for (int icomp=0;icomp<staticthis->nComp;icomp++){
//      for (int iatt=0;iatt<staticthis->nAtt;iatt++){
 //       staticthis->Par[ibin][icomp][iatt][0] = par[index];
 //       staticthis->Par[ibin][icomp][iatt][1] = par[index+1]; 
//        staticthis->thePars->histoPar[ibin][icomp][iatt][0] = par[index];
//        staticthis->thePars->histoPar[ibin][icomp][iatt][1] = par[index+1];
//        index+=2;
//      }
//    }
//  }
  result = (double)staticthis->getTotLnL();
}

void histoCompare::sumSqWrapper(int& ndim, double* gout, double& result, double par[], int flg){
  //index trix
  int index=0;
  for (int ibin=0;ibin<staticthis->nBin;ibin++){
    for (int icomp=0;icomp<staticthis->nComp;icomp++){
      for (int iatt=0;iatt<staticthis->nAtt;iatt++){
        staticthis->Par[ibin][icomp][iatt][0] = par[index];
        staticthis->Par[ibin][icomp][iatt][1] = par[index+1]; 
        index+=2;
      }
    }
  }
  result = (double)staticthis->getTotSumSq();
}

void histoCompare::sumSqDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg){
  //index trix
  int index;
  for (int i=0;i<staticthis->nMCHist;i++){
    index = 2*i;
    staticthis->parDebug[i][0] = par[index];
    staticthis->parDebug[i][1] = par[index+1]; 
  }
  result = (double)staticthis->getTotSumSqDebug();
}

void histoCompare::nDiffDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg){
  //index trix
  int index;
  for (int i=0;i<staticthis->nMCHist;i++){
    index = 2*i;
    staticthis->parDebug[i][0] = par[index];
    staticthis->parDebug[i][1] = par[index+1]; 
  }
  
//  result = (double)staticthis->getTotSumSqDebug();
  result = (double)staticthis->getNDiff();
  
}

void histoCompare::minSumSqDebug(){

  //setup static this
  staticthis = this;

  //setup initial parameters
  int npartot = nMCHist*2; //total parameters in fit
  TFitter* fit = new TFitter(npartot); //fitter obj
  fit->SetFCN(sumSqDebugWrapper);
//  fit->SetFCN(nDiffDebugWrapper);
  for (int ipar=0;ipar<npartot;ipar++){
    if ((ipar%2)==0) fit->SetParameter(ipar,Form("param_%d_smear",ipar),1.0,0.1,0,0);
    else {
      fit->SetParameter(ipar,Form("param_%d_bias",ipar),0.0,0.1,0,0);
    }
  }
  int jpar = 0;
  //fix smear parameters...fit only bias for now
  while (jpar<npartot){
    fit->FixParameter(jpar);
    jpar+=2;
  }
  
  //run the fit
  fit->ExecuteCommand("SIMPLEX",0,0);
  jpar = 0;
  while (jpar<npartot){
    fit->ReleaseParameter(jpar);
    jpar+=2;
  }
  fit->ExecuteCommand("SIMPLEX",0,0);

  
  //get results
  int ihist = 0;
  int index = 0;
  while (index<npartot){
    parDebug[ihist][0] = fit->GetParameter(index);
    parDebug[ihist][1] = fit->GetParameter(index+1);
    cout<<"   $ HISTO: "<<ihist<<" SMEAR: "<<parDebug[ihist][0];
    cout<<" BIAS: "<<parDebug[ihist][1];
    ihist++;
    index+=2;
  } 
  return;
}

void histoCompare::minSumSq(){

  //setup static this
  staticthis = this;
  double parerr = 0.1;
  //setup initial parameters
  int npartot = nBin*nComp*nAtt*2; //total parameters in fit
  int parindex=0;
  TFitter* fit = new TFitter(npartot); //fitter obj
  fit->SetFCN(sumSqWrapper);
  for (int ibin=0;ibin<nBin;ibin++){
    for (int icomp=0;icomp<nComp;icomp++){
      for (int iatt=0;iatt<nAtt;iatt++){
        fit->SetParameter(parindex,Form("param_%d_smear",parindex),Par[ibin][icomp][iatt][0],parerr,0,0);
        if (fixPar[ibin][icomp][iatt][0]) fit->FixParameter(parindex);
        parindex++;
        fit->SetParameter(parindex,Form("param_%d_bias",parindex),Par[ibin][icomp][iatt][1],parerr,0,0);
        if (fixPar[ibin][icomp][iatt][1]) fit->FixParameter(parindex);
        parindex++;
      }
    }
  }  
  //run the fit
  fit->ExecuteCommand("SIMPLEX",0,0);
  //get results
  parindex = 0;
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jcomp=0;jcomp<nComp;jcomp++){
      for (int jatt=0;jatt<nAtt;jatt++){
        Par[jbin][jcomp][jatt][0] = fit->GetParameter(parindex);
        parindex++;
        Par[jbin][jcomp][jatt][1] = fit->GetParameter(parindex);
        parindex++;
      }
    }
  }    
  return;
}

float histoCompare::getTotSumSqDebug(){
  float totsumsq = 0.;
  hMod = smearIt(hMC[0],parDebug[0][0],parDebug[0][1]);
  hTot = (TH1F*)hMod->Clone("hmod");
  for (int ihist=1;ihist<nMCHist;ihist++){
    hTot->Add(smearIt(hMC[ihist],parDebug[ihist][0],parDebug[ihist][1]));
  }
  totsumsq = getSumSq(hTot,hData[0]);
//  for (int i=0;i<nMCHist;i++){
//    hModDebug = smearIt(hTot,parDebug[i][0],parDebug[i][1]);
//    totsumsq+=getSumSq(hModDebug,hData[0]);
//  }
  return totsumsq;
}


float histoCompare::getTotLnL(){
  float totL = 0.;
  for (int isamp=0;isamp<nSamp;isamp++){
    for (int ibin=0;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
         //get modfied MC prediction
      //   if (hMod) hMod->Delete();
        // hMod = (TH1F*)hManager->hMC[isamp][ibin][0][iatt]->Clone("hmodtmp");
   //      hMod = (TH1F*)hManager->getModHistogram(isamp,ibin,0,iatt)->Clone("hmodtmp");
      //   cout<<"getting histogram: "<<isamp<<" "<<ibin<<" "<<iatt<<endl;
         hMod = hManager->getSumHistogramMod(isamp,ibin,iatt);
        // hMod = getModifiedHisto(isamp,ibin,iatt);
     //    hMod->Smooth(5);
      //   hTmp = (TH1F*)hManager->hMC[isamp][ibin][0][iatt]->Clone("htmp");
      //   smearHisto((*hManager->hMC[isamp][ibin][0][iatt]),(*hMod),Par[ibin][0][iatt][0],Par[ibin][0][iatt][1]);    
      //   smearHisto((*hManager->hMC[isamp][ibin][0][iatt]),(*hMod),thePars->histoPars[ibin][0][iatt][0],thePars->histoPars[ibin][0][iatt][1]);    
      //   for (int icomp = 1;icomp<nComp;icomp++){
     //      hTmp->Smooth(1);
    //       smearHisto((*hManager->hMC[isamp][ibin][icomp][iatt]),(*hTmp),Par[ibin][icomp][iatt][0],Par[ibin][icomp][iatt][1]);
    //       hTmp->Smooth(3);
         //  hMod->Add(hManager->getModHistogram(isamp,ibin,icomp,iatt));
        // }
         //add error to total
    //     hMod->Smooth(3);
      //   cout<<"scaling histogram "<<endl;
         //hMod->Scale(Norm);
        // hMod->Smooth(10);
     //    convolveThisHisto(*hMod,hMod->GetBinWidth(2)*0.7,0.);
     //    cout<<"rebin histogram with factor: "<<rebinFactor<<endl;
       //  hMod->Rebin(rebinFactor);
         totL+=getLnL(hMod,hManager->hData[isamp][ibin][iatt]);
      }
    }
  }
  
   //add flux/xsec component
  float pull;
  for (int isys=0;isys<thePars->nSysPars;isys++){
    pull = thePars->sysPar[isys]-1.;
    pull/=thePars->sysParUnc[isys];
    totL+=(0.5)*(pull*pull);
  }
  return totL;
}



float histoCompare::getTotSumSq(){
  float totsumsq = 0.;
  for (int isamp=0;isamp<nSamp;isamp++){
    for (int ibin=0;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
         //get modfied MC prediction
         hMod = smearIt(hManager->hMC[isamp][ibin][0][iatt],Par[ibin][0][iatt][0],Par[ibin][0][iatt][1]);      
         for (int icomp = 1;icomp<nComp;icomp++){
           hMod->Add(smearIt(hManager->hMC[isamp][ibin][icomp][iatt],Par[ibin][icomp][iatt][0],Par[ibin][icomp][iatt][1]));
         }
         //add error to total
        // hMod->Scale(Norm);
         totsumsq+=getSumSq(hMod,hManager->hData[isamp][ibin][iatt]);
      }
    }
  }
  return totsumsq;
}


float histoCompare::getSumSq(TH1F* h1, TH1F* h2){
  float sumsq = 0.;
  float diff;
  for (int ibin=10;ibin<=(h1->GetNbinsX()-10);ibin++){
    diff = h1->GetBinContent(ibin)-h2->GetBinContent(ibin);
    sumsq += (diff*diff);
  }
  return sumsq;
}

float histoCompare::getLnL(TH1F* h1, TH1F* h2){
  float lnL = 0.;
  float diff;
  float term;
  float c1;
  float c2;
  for (int ibin=3;ibin<=(h1->GetNbinsX()-3);ibin++){
    c1 = h1->GetBinContent(ibin);
    c2 = h2->GetBinContent(ibin);
    diff = c1-c2;
   // if (c1==0) return (c1-c2)*(c1-c2);
    if (c2==0){
      term = 0.;
    }
    else if(c1==0.){
      term = -1.*diff;
    }
    else{ 
      term = c2*TMath::Log(c2/c1);
    }
    lnL += (diff+term);
  }
  return lnL;
}

void histoCompare::addHistogram(TH1F* h,int dataflg){
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

void histoCompare::readFromFile(const char* filerootname,int nsamp, int nbin, int ncomp, int natt){
  nSamp = nsamp;
  nBin  = nbin;
  nComp = ncomp; 
  nAtt  = natt;
  hManager = new histoManager(filerootname,nsamp,nbin,ncomp,natt);
  float ndataevents=0;
  float nmcevents=0;
  //setup initaal parametres
  /*
  for (int icomp=0;icomp<nComp;icomp++){
    for (int ibin=0;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
        Par[ibin][icomp][iatt][0]=1.;  //smear parameter
        Par[ibin][icomp][iatt][1]=0.;  //bias parameter
        bestPar[ibin][icomp][iatt][0]=1.;  //smear parameter
        bestPar[ibin][icomp][iatt][1]=0.;  //bias parameter
      }
    }
  } 
  */
  float events;
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
  Norm = ndataevents/nmcevents;
  
  return;
}

void histoCompare::setupPars(int nsyspars){
  thePars = new atmFitPars(nSamp,nBin,nComp,nAtt,"tn186");
 // thePars->setNorm(Norm);
  hManager->setFitPars(thePars);
  return;
}


histoCompare::histoCompare(const char* thename){
  nameTag = thename;
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
  return;
}
#endif
