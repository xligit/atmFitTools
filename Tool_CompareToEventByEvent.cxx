#include "histoManager.cxx"
#include "hSplines.cxx"
#include "atmFitPars.cxx"

using namespace std;


class compareToEventByEvent{
  public:
  compareToEventByEvent(atmFitPars* fitpars, TTree* tr, const char* hfilename, const char* sfilename);
  histoManager* hManager;
  TH1F* htmp1;
  TH1F* htmp2;
  TH1F* htmp3;
  TTree* mcTree;
  atmFitPars* thePars;
  fQreader* mcEvt;
  float getEvtWeight(int isys);
  float att[1000];
  void fillAttributes();
  void comparePrediction(int isamp,int ibin, int icomp, int iatt, int isys);
};

void compareToEventByEvent::comparePrediction(int isamp,int ibin,int icomp,int iatt,int isys){
  //clone histogram
  htmp1 = (TH1F*)hManager->getHistogram(isamp,ibin,icomp,iatt)->Clone("tmpclone");
  htmp1->Reset();
  htmp2 = hManager->getModHistogram(isamp,ibin,icomp,iatt);
  htmp3 = hManager->getHistogram(isamp,ibin,icomp,iatt);
  
  //fill histogram from tree
  for (int iev=0;iev<mcTree->GetEntries();iev++){
    mcTree->GetEvent(iev);
    fillAttributes();
    float weight = getEvtWeight(isys);
    if ((mcEvt->nbin==ibin)&&(mcEvt->nsample==isamp)&&(mcEvt->ncomponent==icomp)){
      htmp1->Fill(att[iatt],weight);
    }
  }
  htmp1->SetLineColor(kBlack);
  htmp2->SetLineColor(kRed);
  htmp3->SetLineColor(kBlue);
  htmp1->Draw("");
  htmp1->Draw("samee");
  htmp2->Draw("samee");
  htmp2->Draw("sameh");
  htmp3->Draw("samee");
  htmp3->Draw("sameh");
  return;
}

void compareToEventByEvent::fillAttributes(){
  att[0] = mcEvt->fq1rnll[0][2]-mcEvt->fq1rnll[0][1];
  att[1] = mcEvt->fq1rnll[1][2]-mcEvt->fq1rnll[1][1];
  return;
}

float compareToEventByEvent::getEvtWeight(int ipar){
//  float ww = 1.;
  float ww = mcEvt->evtweight;
  int absmode = TMath::Abs(mcEvt->mode);
  float Enu     = mcEvt->pmomv[0];
  int  nutype  = TMath::Abs(mcEvt->ipnu[0]);
   
  //CCQE norm bin1 
 // if (ipar==0){
    if ((absmode==1)&&(Enu<200.)) ww*=thePars->sysPar[0];
 // }
  //CCQE norm bin2 
//  if (ipar==1){
    if ((absmode==1)&&(Enu>200.)&&(Enu<400.)) ww*=thePars->sysPar[1];
 // }
  //CCQE norm bin3 
 // if (ipar==2){
    if ((absmode==1)&&(Enu>400.)&&(Enu<800.)) ww*=thePars->sysPar[2];
//  }
  //CCQE norm bin4 
 // if (ipar==3){
    if ((absmode==1)&&(Enu>800.)) ww*=thePars->sysPar[3];
 // }
  //SubGevFlux
 // if (ipar==4){
    if (Enu<1000.) ww*=thePars->sysPar[4];
 // }
  //MultiGeVFlux
//  if (ipar==5){
    if (Enu>1000.) ww*=thePars->sysPar[5];
//  }
  //CCnQE
//  if (ipar==6){
    if ((absmode>1)&&(absmode<30)) ww*=thePars->sysPar[6];
//  }
  //NC
//  if (ipar==7){
    if (absmode>=30) ww*=thePars->sysPar[7];
//  }
  //mu2e ratio
//  if (ipar==8){
    if (nutype==14) ww*=thePars->sysPar[8];
//  }

  if (ww<0.) ww = 0.;
  return ww;
}

compareToEventByEvent::compareToEventByEvent(atmFitPars* thepars, TTree* tr,const char* hfilename, const char* sfilename){
  //setup mc tree
  mcTree = tr;
  mcEvt = new fQreader(mcTree);
  //setup fit pars
  thePars = thepars;
  //creaate manager
  hManager = new histoManager(hfilename,thePars->nSamples,thePars->nBins,thePars->nComponents,thePars->nAttributes);
  //read in spline info
  hManager->setFitPars(thePars);
  hManager->readSplinesFromFile(sfilename,thePars->nSysPars); 
}


