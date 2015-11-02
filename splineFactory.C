#include "hSplines.C"
#include "shared.h"
//#include "fQreader.C"
#include "histoManager.C"

#define NSYSTMAX 10
#define NHBINSMAX 300
#define NPTSMAX   5

using namespace std;

//class for creating splines
class splineFactory{
  public:
  //constructor
  splineFactory(int nsamp, int nbin, int ncomp, int natt, int nsyst, const char* name);
  splineFactory(){;}
  

  //internal vars
  TString nameTag;
  TTree* mcTree;
  fQreader* mcEvt;
  histoManager* hManager; //manages all default histograms

  TH1F* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; //array for modified
  TFile *fout;
  int nSamp;
  int nBin;
  int nComp;
  int nAtt;
  int nSyst;
  float sysPar[NSYSTMAX];
  float sysUnc[NSYSTMAX];
  float attribute[NATTMAX];
  float eventWeight;

  //for output tree
  TTree* splineTree;
  int nbin;
  int ncomponent;
  int nattribute;
  int nsample;
  int nsystpar;
  int npoints;
  int nhistobins;
  float systParValues[NPTSMAX];
  float binWeight[NPTSMAX][NHBINSMAX];
  //methods
  float getEvtWeight(int ipar); //returns event weight after applying syst. par. 
  void fillHistograms(int ipt,int isyst);
  void  makeManagerFromFile(const char* fname);
  void fillLeaves(int nsamp,int nbin,int ncomp,int natt,int isyst);  
  void setMCTree(TTree* tr);
  //build the splines
  void buildTheSplines();

  //debugging
  void debugtest();

//  private:
  void fillAttributes();
  void setupHistos();
  void setupSystPars();
  void incrementSystPars(float nsig);
};


void splineFactory::setMCTree(TTree* tr){
  mcTree = tr;
  mcEvt = new fQreader(mcTree);
  return;
}


void splineFactory::debugtest(){
  //create new factory
  splineFactory* sfactory = new splineFactory(3,3,7,1,1,"debugtest");
  //set histogram manager and create template histograms
  sfactory->makeManagerFromFile("factoryOut_factorytest.root");
  

  return;
};

void splineFactory::fillLeaves(int isamp,int ibin,int icomp,int iatt,int isyst){
   nsample = isamp;
   nhistobins = hMC[isamp][ibin][icomp][iatt][0]->GetNbinsX();
   npoints = NPTSMAX;
   for (int ipt=0;ipt<NPTSMAX;ipt++){
     for (int jhistobin=0;jhistobin<=nhistobins;jhistobin++){
       if (hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin)==0){
         binWeight[ipt][jhistobin] = 1.;
       }
       else{
         binWeight[ipt][jhistobin] =
           (float)hMC[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin)/
           (float)hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin);
       }
       cout<<"value: "<<ipt<<","<<jbin<<" -> "<<binWeight[ipt][jhistobin]<<endl;
     }
   }
   return;
}



void splineFactory::buildTheSplines(){
  // file setup
  TString fname = "splineOut_";
  fname.Append(nameTag.Data());
  fname.Append(".root");
  fout = new TFile(fname.Data(),"recreate");

  //setup the output tree
  splineTree = new TTree("splinePars","spinePars");
  splineTree->Branch("nbin",&nbin,"nbin/I");
  splineTree->Branch("nhistobins",&nhistobins,"nhistobins/I");
  splineTree->Branch("ncomponent",&ncomponent,"ncomponent/I");
  splineTree->Branch("nattribute",&nattribute,"nattribute/I");
  splineTree->Branch("nsample",&nsample,"nsample/I");
  splineTree->Branch("nsystpar",&nsystpar,"nsystpar/I");
  splineTree->Branch("nsyspartot",&nSyst,"nsyspartot/I");
  splineTree->Branch("npoints",&npoints,"npoints/I");
  splineTree->Branch("systParValues",systParValues,Form("systParValues[%d]/F",NPTSMAX));
  splineTree->Branch("binWeight",binWeight,Form("binWeight[%d][%d]/F",NPTSMAX,NHBINSMAX));

  

  //setup systematic deviations
  float sigvals[5] = {-5.,-2.,0.,2.,5.};
  cout<<"creating spines"<<endl; 
  for (int isyst=0;isyst<nSyst;isyst++){
    for (int ipt=0;ipt<NPTSMAX;ipt++){
      //fix the systematic parameters for this point
      incrementSystPars(sigvals[ipt]);
      //fill array of systematic parameter values
      systParValues[ipt]=sysPar[isyst];
      //loop over all events and fill re-weighted histograms
      cout<<"filling histograms for syspar: "<<isyst<<" point: "<<ipt<<endl;
      fillHistograms(ipt,isyst);
    }
    //loop over all histograms and get spline information
    //and write to tree
    for (int kbin=0;kbin<nBin;kbin++){
      for (int kcomp=0;kcomp<nComp;kcomp++){
        for (int ksamp=0;ksamp<nSamp;ksamp++){
          for (int katt=0;katt<nAtt;katt++){
            fillLeaves(ksamp,kbin,kcomp,katt,isyst);
            splineTree->Fill();
          }
        }
      }
    }
  }
  splineTree->Write();
  fout->Close();
  return;
}

void splineFactory::incrementSystPars(float nsig){

  //reset initial parameters
  setupSystPars();

  //change systematic parameters
  for (int isyst=0;isyst<nSyst;isyst++){
    sysPar[isyst] += sysUnc[isyst]*nsig;
  }
  
  return;
}

void splineFactory::setupSystPars(){
  nSyst=1;
  sysPar[0] = 1.0;
  sysUnc[0] = 0.1;  //one sigma uncertainty
  return;
}

void splineFactory::setupHistos(){
  cout<<"called setupHistos()"<<endl;
  TH1F* htemplate;
  TString hname;
  for (int ipt=0;ipt<NPTSMAX;ipt++){
    for (int isamp=0;isamp<nSamp;isamp++){
      for (int ibin=0;ibin<nBin;ibin++){
        for (int icomp=0;icomp<nComp;icomp++){
          for (int iatt=0;iatt<nAtt;iatt++){
            htemplate = hManager->getHistogram(isamp,ibin,icomp,iatt);
            hname = htemplate->GetName();
            hname.Append(Form("_%d%d%d%d%d",ipt,isamp,ibin,icomp,iatt));
            cout<<"setting histogram template: "<<hname.Data()<<endl;
            hMC[isamp][ibin][icomp][iatt][ipt]=(TH1F*)htemplate->Clone(hname.Data());
            hMC[isamp][ibin][icomp][iatt][ipt]->Reset();
          }
        }
      }
    }
  }
  return;
}

void splineFactory::makeManagerFromFile(const char* fname){
  hManager = new histoManager(fname,nSamp,nBin,nComp,nAtt);
  //setupHistos();
  return;
}

void splineFactory::fillAttributes(){
  attribute[0] = mcEvt->fq1rnll[0][2]-mcEvt->fq1rnll[0][1];
  attribute[1] = mcEvt->fq1rnll[1][2]-mcEvt->fq1rnll[1][1];
  return;
}

void splineFactory::fillHistograms(int ipt, int isyst){
  //fills histograms after applying systematic error parameter
  for (int iev=0;iev<mcTree->GetEntries();iev++){
     mcTree->GetEvent(iev);
     fillAttributes();
     getEvtWeight(isyst);
     for (int jatt=0;jatt<nAtt;jatt++){
      hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][jatt][ipt]->Fill(attribute[jatt],eventWeight);
    }
  }

  return;
};

float splineFactory::getEvtWeight(int ipar){
  float ww = 1.;
  if (ipar==0){
    if (TMath::Abs(mcEvt->mode)<30.) ww*=sysPar[0];
  }
  eventWeight = ww;
  return ww;
};

splineFactory::splineFactory(int isamp, int ibin, int icomp, int iatt, int isyst, const char* name){
  nameTag = name;
  nSamp = isamp;
  nBin  = ibin;
  nComp = icomp;
  nAtt  = iatt;
  nSyst = isyst;
}

