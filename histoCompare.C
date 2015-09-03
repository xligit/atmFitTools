#include "histoManager.C"
#include "TMath.h"
#include "TRandom2.h"
#include "smearHisto.C"
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
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp);
  histoManager* hManager;
 // TH1F* hMod[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];

  //tools for adding histogram directly..for debugging
  void addHistogram(TH1F* h,int dataflg);
  TH1F* hData[10];
  TH1F* hMC[10];
  TH1F* hModDebug;
  TH1F* hMod;
  void showMod(int imchist);
  TH1F* hTot;
  int nDataHist;
  int nMCHist;
  float parDebug[10][2];
//  float getLDebug(); 
  
  //likelihood evaluateions
  float getSumSq(TH1F* h1, TH1F* h2);
  float getTotSumSqDebug();
  static void sumSqDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void minSumSqDebug();
  void drawResult(int ihist);

  //staticthis for fits
  static histoCompare* staticthis;
};

histoCompare* histoCompare::staticthis;

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

void histoCompare::minSumSqDebug(){

  //setup static this
  staticthis = this;

  //setup initial parameters
  int npartot = nMCHist*2; //total parameters in fit
  TFitter* fit = new TFitter(npartot); //fitter obj
  fit->SetFCN(sumSqDebugWrapper);
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

float histoCompare::getSumSq(TH1F* h1, TH1F* h2){
  float sumsq = 0.;
  float diff;
  for (int ibin=5;ibin<=(h1->GetNbinsX()-5);ibin++){
    diff = h1->GetBinContent(ibin)-h2->GetBinContent(ibin);
    sumsq += (diff*diff);
  }
  return sumsq;
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

void histoCompare::readFromFile(const char* filerootname,int isamp, int ibin, int icomp){
  nSamp = isamp;
  nBin  = ibin;
  nComp = icomp; 
  TString hmname = "hManagerFor_";
  hmname.Append("nametag.data");
  hManager = new histoManager(isamp,ibin,icomp,hmname.Data());
  hManager->readFromFile(filerootname);
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
