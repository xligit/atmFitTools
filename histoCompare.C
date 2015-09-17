#include "histoManager.C"
#include "TMath.h"
#include "TRandom2.h"
#include "histoTransforms.C"
#include "TFitter.h"

using namespace std;


//class to compare histograms and evaluate likelihoods
class histoCompare{
  public:

  //constructors//
  histoCompare(const char* thename);  //standard constructor

  //internal variables
  TString nameTag;  //name associated with this instance
  int nSamp;  //number of samples
  int nBin;  //number of bins
  int nComp;  //number of  components
  int nAtt;  //nummboer of attributes
  //tools for histogram manager management
  //created histo manager from file
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp, int natt);
  histoManager* hManager;
  float Norm;
  float Par[NBINMAX][NCOMPMAX][NATTMAX][2];
  float fixPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  float bestPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  float errPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString parName[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString binName[NBINMAX];
  TString compName[NCOMPMAX];
  TString attName[NCOMPMAX];
  void setBinName(int ibin, const char* name){binName[ibin]=name;}
  void setCompName(int icomp, const char* name){compName[icomp]=name;}
  void setAttName(int iatt, const char* name){attName[iatt]=name;}
  void showFitHisto(int isamp,int ibin,int icomp,int iatt);
  void showFitEffect(int isamp,int ibin,int icomp,int iatt);
  void showFitResult(int isamp,int ibin,int iatt);
  void showFitPars(int ibin,int iatt,int imod);
 // TH1F* hMod[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];

  //tools for adding histogram directly..for debugging
  void addHistogram(TH1F* h,int dataflg);
  TH1F* hData[10];
  TH1F* hMC[10];
  TH1F* hModDebug;
  TH1F* hMod;
  TH1F* hPar;
  TH1F* showSmear(TH1F* h, float smear, float bias);
  void showMod(int imchist);
  TH1F* hTot;
  int nDataHist;
  int nMCHist;
  float parDebug[10][2];
//  float getLDebug(); 
  
  //likelihood evaluateions
  float getSumSq(TH1F* h1, TH1F* h2);
  float getLnL(TH1F* h1, TH1F* h2);
  float getNDiff();
  float getTotSumSq();
  float getTotLnL();
  static void sumSqWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  //for debuggint and play
  float getTotSumSqDebug();
  static void sumSqDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void nDiffDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void minSumSqDebug();
  void minSumSq();
  void sumSqPrefit();
  void LnLFit();
  void LnLPreFit();
  void drawResult(int ihist);

  //staticthis for fits
  static histoCompare* staticthis;
};

histoCompare* histoCompare::staticthis;

void histoCompare::showFitPars(int ibin,int iatt,int imod){
  int nbinstot = nComp;
  hPar  = new TH1F("hpar","hpar",nbinstot,0,nbinstot);
  for (int icomp=0;icomp<nComp;icomp++){
    hPar->SetBinContent(icomp+1,bestPar[ibin][icomp][iatt][imod]);
    hPar->GetXaxis()->SetBinLabel(icomp+1,parName[ibin][icomp][iatt][imod].Data());
  }
  hPar->Draw();
  return;
}


TH1F* histoCompare::showSmear(TH1F* h, float smear, float bias){
  TH1F* hh = smearIt(h,smear,bias);
  return hh;
}


void histoCompare::showFitResult(int isamp,int ibin,int iatt){
  //sum up mc components
  float smear = bestPar[ibin][0][iatt][0];
  float bias = bestPar[ibin][0][iatt][1];
  hMod = smearIt(hManager->hMC[isamp][ibin][0][iatt],smear,bias);
  hTot = (TH1F*)hManager->hMC[isamp][ibin][0][iatt]->Clone("htot");
  for (int jcomp=1;jcomp<nComp;jcomp++){
    hTot->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
    smear = bestPar[ibin][jcomp][iatt][0];
    bias  = bestPar[ibin][jcomp][iatt][1];
    hMod->Add(smearIt(hManager->hMC[isamp][ibin][jcomp][iatt],smear,bias));
  }
  hTot->SetLineColor(kRed);
  hTot->Scale(Norm);
  hMod->SetLineColor(kBlue);
  hMod->Scale(Norm);
  hTot->Draw();
  hMod->Draw("same");
  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->Draw("samee");
  return;
}

void histoCompare::showFitEffect(int isamp,int ibin,int icomp,int iatt){
  //sum up mc components
  float smear = bestPar[ibin][icomp][iatt][0];
  float bias  = bestPar[ibin][icomp][iatt][1];
  cout<<"SMEAR: "<<smear<<endl;
  cout<<"BIAS:  "<<bias<<endl;
  hMod = smearIt(hManager->hMC[isamp][ibin][icomp][iatt],smear,bias);
  hTot = (TH1F*)hManager->hMC[isamp][ibin][icomp][iatt]->Clone("htot");
  for (int jcomp=0;jcomp<nComp;jcomp++){
    if (jcomp!=icomp){
      hMod->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
      hTot->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
    }
  }
  hTot->SetLineColor(kRed);
  hTot->Scale(Norm);
  hMod->SetLineColor(kBlue);
  hMod->Scale(Norm);
  hTot->Draw();
  hMod->Draw("same");
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
  hMod->Draw("same");
  return;
}

void histoCompare::LnLPreFit(){
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
  int parindextmp;
  //parameter name container
  TString parnametmp;  

  float  parinit; //container to temporarily store initial values 

  //total number of parameters to be fit
  int npars = nBin*nComp*nAtt*2;

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
    double pp = 1;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //set parameters to inital values
  parindex = 0;
  for (int kbin=0;kbin<nBin;kbin++){
    for (int katt=0;katt<nAtt;katt++){
      for (int kcomp=0;kcomp<nComp;kcomp++){
        bestPar[kbin][kcomp][katt][0] = Par[kbin][kcomp][katt][0];
        fit->SetParameter(parindex,parName[kbin][kcomp][katt][0].Data(),Par[kbin][kcomp][katt][0],parerr,0,0);
        parindex++;
        bestPar[kbin][kcomp][katt][1] = Par[kbin][kcomp][katt][1];
        fit->SetParameter(parindex,parName[kbin][kcomp][katt][1].Data(),Par[kbin][kcomp][katt][1],parerr,0,0);
        parindex++;
      }
    }
  }

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }

  parindex = 0;
  //run individual fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
        //fix all other variables
        //release parameter to be fit
        parindex++;
        fit->ReleaseParameter(parindex);
        fit->ExecuteCommand("SIMPLEX",0,0);
        fit->FixParameter(parindex);
        bestPar[jbin][jcomp][jatt][1] = fit->GetParameter(parindex); 
        parindex++;
      }
    }
  }

  //print final results
  
  for (int pbin=0;pbin<nBin;pbin++){
    for (int patt=0;patt<nAtt;patt++){
      for (int pcomp=0;pcomp<nComp;pcomp++){
        for (int pmod=0;pmod<2;pmod++){
          cout<<"  PAR "<<parName[pbin][pcomp][patt][pmod].Data()<<" FIT RESULT: ";
          cout<<bestPar[pbin][pcomp][patt][pmod]<<endl;
        }
      }
    }
  }

  cout<<"$$$$$$$$$$$$$ END LNL PRE FIT $$$$$$$$$$$$$"<<endl;

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
  int npars = nBin*nComp*nAtt*2;

  cout<<"$$$$$$$$$$$$$$$$$ LNL FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //run the prefit
  LnLPreFit();

  //setup the fitter!
  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = 1;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //set parameters to inital values (from prefit)
  parindex = 0;
  for (int kbin=0;kbin<nBin;kbin++){
    for (int katt=0;katt<nAtt;katt++){
      for (int kcomp=0;kcomp<nComp;kcomp++){
        fit->SetParameter(parindex,parName[kbin][kcomp][katt][0].Data(),bestPar[kbin][kcomp][katt][0],parerr,0,0);
        parindex++;
        fit->SetParameter(parindex,parName[kbin][kcomp][katt][1].Data(),bestPar[kbin][kcomp][katt][1],parerr,0,0);
        parindex++;
      }
    }
  }
  
  parindex = 0;
  int parindextmp = 0;
  //do individual fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      //start of fit block
      //fix all parameters
      for (int jpar=0;jpar<npars;jpar++){
        fit->FixParameter(jpar);
      }
      parindextmp = parindex;
      //release bias parameters
      for (int jcomp=0;jcomp<nComp;jcomp++){
        //fill initial value to restore parameter array later
        bestPar[jbin][jcomp][jatt][1] = Par[jbin][jcomp][jatt][1];
        //fit a single parameter
        fit->ReleaseParameter(parindex+1);   
        parindex+=2;
      }

      //run fit
      fit->ExecuteCommand("SIMPLEX",0,0); //run the fit for bias parameters
      //release smear parameters
      parindex = parindextmp;
      for (int jcomp=0;jcomp<nComp;jcomp++){
        //fill initial value to restore parameter array later
        bestPar[jbin][jcomp][jatt][0] = Par[jbin][jcomp][jatt][0];
        //fit a single parameter
        fit->ReleaseParameter(parindex);   
        parindex+=2;
      }
      fit->ExecuteCommand("SIMPLEX",0,0); //run the fit for ALL parameters

      parindex = parindextmp;
      //get results
      for (int lcomp=0;lcomp<nComp;lcomp++){
        for (int lmod=0;lmod<2;lmod++){
          Par[jbin][lcomp][jatt][lmod] = bestPar[jbin][lcomp][jatt][lmod]; 
          bestPar[jbin][lcomp][jatt][lmod] = fit->GetParameter(parindex); 
          errPar[jbin][lcomp][jatt][lmod] = fit->GetParError(parindex);
          parindex++;
        }
      }
      //end of fit block
    }
  }

  //print final results
  cout<<"  ----------------------------------------  "<<endl;
  for (int pbin=0;pbin<nBin;pbin++){
    for (int patt=0;patt<nAtt;patt++){
      for (int pcomp=0;pcomp<nComp;pcomp++){
        for (int pmod=0;pmod<2;pmod++){
          cout<<"  PAR "<<parName[pbin][pcomp][patt][pmod].Data()<<" FIT RESULT: ";
          cout<<Par[pbin][pcomp][patt][pmod]<<" -> "<<bestPar[pbin][pcomp][patt][pmod]<<endl;
    //      cout<<" +/- "<<errPar[pbin][pcomp][patt][pmod]<<endl;
        }
      }
    }
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

void histoCompare::lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg){
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
         hMod = smearIt(hManager->hMC[isamp][ibin][0][iatt],Par[ibin][0][iatt][0],Par[ibin][0][iatt][1]);      
         for (int icomp = 1;icomp<nComp;icomp++){
           hMod->Add(smearIt(hManager->hMC[isamp][ibin][icomp][iatt],Par[ibin][icomp][iatt][0],Par[ibin][icomp][iatt][1]));
         }
         //add error to total
         hMod->Scale(Norm);
         totL+=getLnL(hMod,hManager->hData[isamp][ibin][iatt]);
      }
    }
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
         hMod->Scale(Norm);
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
  for (int ibin=10;ibin<=(h1->GetNbinsX()-10);ibin++){
    c1 = h1->GetBinContent(ibin);
    c2 = h2->GetBinContent(ibin);
    diff = c1-c2;
    term = c2*TMath::Log(c2/c1);
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
  //count total events
  for (int jsamp=0;jsamp<nSamp;jsamp++){
    for (int jbin=0;jbin<nBin;jbin++){
      for (int jatt=0;jatt<nAtt;jatt++){
        ndataevents+=(float)hManager->hData[jsamp][jbin][jatt]->Integral();
        for (int jcomp=0;jcomp<nComp;jcomp++){
          nmcevents+=(float)hManager->hMC[jsamp][jbin][jcomp][jatt]->Integral();
        }
      }
    }
  }
  Norm = ndataevents/nmcevents;
  return;
}

histoCompare::histoCompare(const char* thename){
  nameTag = thename;
  nMCHist=0;
  nDataHist=0;
  cout<<"created comparison object: "<<nameTag.Data()<<endl;
  //setup initial debug
  int jhist = 0;
  while (jhist<10){
      parDebug[jhist][0] = 1.0;
      parDebug[jhist][1] = 0.0;
      jhist++;
  }
  return;
}
