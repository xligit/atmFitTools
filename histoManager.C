#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "fQreader.C"
#include <iostream>

#define NSAMPMAX 5
#define NCOMPMAX 20
#define NATTMAX 20
#define NBINMAX 10

class histoManager{
  public:
  histoManager(int nsampl,int nbins,int ncomp,const char* name=""); //constructor
  histoManager(const char* rootname);
  TTree* dataTree;
  TTree* mcTree; 
  fQreader* fqData;
  fQreader* fqMC;
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TFile*   fin; //input file of histos to be read
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  TH1F* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TH1F* hData[NSAMPMAX][NBINMAX][NATTMAX];  //array of all Data histograms
  float att[NATTMAX]; //array of all attribute values
  TString attNames[NATTMAX];  //array of attribute names
  TString attType[NATTMAX];  //array of attribute type codes
  void init();  //initialize after attributes have been set
  void addAttribute(int iatt);  //add an attribute (fiTQun variable) to list
  
  //histogram maker
  TH1F* getHistogram(int iatt,const char* thename); 
  
  //fill attribute variables
  void fillAttributesData();
  void fillAttributesMC();
  void fillHistos();

  //setters
  void setDataTree(TTree* tr);
  void setDataTree(TChain* ch);
  void setMCTree(TTree* tr);
  void setMCTree(TChain* ch);
 
  //file management
  void saveToFile();
  void readFromFile(const char* rootename);
};

histoManager::histoManager(const char* rootname){
  readFromFile(rootname);
  return;
}

void histoManager::readFromFile(const char* rootname){
  TString filename = rootname;
  filename.Append(".root");
  fin = new TFile(filename.Data());
  TString hname;
  //setup data histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int iatt=0;iatt<nAttributes;iatt++){
         hname = "hdata_";
         hname.Append(Form("samp%d_bin%d_att%d",isamp,ibin,iatt));
         cout<<"Getting histogram: "<<hname.Data()<<endl;
         hData[isamp][ibin][iatt] = (TH1F*)fin->Get(hname.Data());
      }
    }
  }
  //setup mc histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int icomp=0;icomp<nComponents;icomp++){
        for (int iatt=0;iatt<nAttributes;iatt++){
           hname = "hmc_";
           hname.Append(Form("samp%d_bin%d_comp%d_att%d",isamp,ibin,icomp,iatt));
           cout<<"Getting histogram: "<<hname.Data()<<endl;
           hMC[isamp][ibin][icomp][iatt] = (TH1F*)fin->Get(hname.Data());
        }
      }
    }
  } 
  return;
}

void histoManager::setMCTree(TChain* ch){
  mcTree=(TTree*)ch;
  fqMC = new fQreader(mcTree);
  return;
}

void histoManager::setMCTree(TTree* tr){
  mcTree=tr;
  fqMC = new fQreader(mcTree);
  return;
}

void histoManager::setDataTree(TChain* ch){
  dataTree=(TTree*)ch;
  fqData = new fQreader(dataTree);
  return;
}

void histoManager::setDataTree(TTree* tr){
  dataTree=tr;
  fqData = new fQreader(dataTree);
  return;
}

void histoManager::saveToFile(){
  fout->Write();
  return;
}

void histoManager::fillHistos(){
  int nevdata = dataTree->GetEntries();
  int nevmc   = mcTree->GetEntries();
  //fill data histos
  for (int i=0;i<nevdata;i++){
    if ((i%100)==0) cout<<"getting data event: "<<i<<endl;
    dataTree->GetEntry(i);
    fillAttributesData();
    for (int iatt=0;iatt<nAttributes;iatt++){
      hData[fqData->nsample][fqData->nbin][iatt]->Fill(att[iatt]);
   }
  }
  for (int j=0;j<nevmc;j++){
    if ((j%100)==0) cout<<"getting mc event: "<<j<<endl;
    mcTree->GetEntry(j);
    fillAttributesMC();
    for (int jatt=0;jatt<nAttributes;jatt++){
      hMC[fqMC->nsample][fqMC->nbin][fqMC->ncomponent][jatt]->Fill(att[jatt]);
    }
  } 
  return;
}

void histoManager::fillAttributesData(){
  att[0] = fqData->fq1rnll[0][2]-fqData->fq1rnll[0][1];
  att[1] = fqData->fq1rnll[1][2]-fqData->fq1rnll[1][1];
  return;
}

void histoManager::fillAttributesMC(){
  att[0] = fqMC->fq1rnll[0][2]-fqMC->fq1rnll[0][1];
  att[1] = fqMC->fq1rnll[1][2]-fqMC->fq1rnll[1][1];
  return;
}


TH1F* histoManager::getHistogram(int iatt, const char* thename){
  TH1F* hnew;
  int nBinsNllEMu = 100;
  if (iatt==0){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-10000,2000);
  }
  if (iatt==1){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-10000,2000);
  }
  return hnew;
}

void histoManager::init(){
  //call ONLY after all attributes are specified
  TString fname = nameTag.Data();
  fname.Append(".root");
  fout = new TFile(fname.Data(),"RECREATE");
  //setup attribute names
  attNames[0]="fq1rnll_emu_subev0";
  attNames[1]="fq1rnll_emu_subev1";
  TString hname;
  //setup data histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int iatt=0;iatt<nAttributes;iatt++){
         hname = "hdata_";
         hname.Append(Form("samp%d_bin%d_att%d",isamp,ibin,iatt));
         cout<<"Making histogram: "<<hname.Data()<<endl;
         hData[isamp][ibin][iatt] = getHistogram(iatt,hname.Data()); //make taylored histo
      }
    }
  }
  //setup mc histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int icomp=0;icomp<nComponents;icomp++){
        for (int iatt=0;iatt<nAttributes;iatt++){
           hname = "hmc_";
           hname.Append(Form("samp%d_bin%d_comp%d_att%d",isamp,ibin,icomp,iatt));
           cout<<"Making histogram: "<<hname.Data()<<endl;
           hMC[isamp][ibin][icomp][iatt] = getHistogram(iatt,hname.Data()); //make taylored histo
        }
      }
    }
  }
  return;
}


void histoManager::addAttribute(int iatt){
  attType[nAttributes]=iatt;
  nAttributes++;
  return;
}

histoManager::histoManager(int nsampl,int nbins,int ncomp,const char* name){
  nameTag = "hManager";
  nameTag.Append(name);
  nSamples = nsampl;
  nComponents = ncomp;
  nAttributes = 0;
  nBins = nbins;
  return;
}


