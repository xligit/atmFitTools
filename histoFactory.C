#include "histoFactory.h"

//constructor to re-created a histogram factory from a file
histoFactory::histoFactory(int nsampl,int nbins,int ncomp,const char* name){
  nameTag = name;
  nameTag.Append("_hFactoryOutput");
  nSamples = nsampl;
  nComponents = ncomp;
  nAttributes = 0;
  nBins = nbins;
  //create an empty histogram manager
  hManager = new histoManager(nsampl,nbins,ncomp,name); 
  return;
}

void histoFactory::addAttribute(int iatt){
  //adds an attribute to the attribute list
  //the attribute type is a code for the fiTQun output being used
  //this code is used to determine histogram binning and names

  //list of codes:
  //0 - e/mu likelihood ratio for subev 1
  //1 - e/mu likelihood ratio for subev 2
  //

  attType[nAttributes]=iatt; 
  nAttributes++;
  return;
}


void histoFactory::init(){
  //call ONLY after all attributes are specified
  //setup name for output file
  TString fname = nameTag.Data();
  fname.Append(".root");
  outputFileName=fname.Data();
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
         hManager->setHistogram(isamp,ibin,0,iatt,1,getHistogramData(iatt,hname.Data()));
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
           hManager->setHistogram(isamp,ibin,icomp,iatt,0,getHistogramData(iatt,hname.Data()));
        }
      }
    }
  }
  return;
}


//use this function to build histograms
TH1F* histoFactory::getHistogramData(int iatt, const char* thename){
  //histogram building
  TH1F* hnew;
  int nBinsNllEMu = 50;
  if (iatt==0){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-3000,6000);
  }
  if (iatt==1){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-3000,6000);
  }
  return hnew;
}



TH1F* histoFactory::getHistogram(int iatt, const char* thename){
  TH1F* hnew;
  int nBinsNllEMu = 50;
  int nBinsNllEMuData = 50;
  if (iatt==0){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-3000,6000);
  }
  if (iatt==1){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-3000,6000);
  }
  return hnew;
}

//calculates attributes from the raw fiTQun output
void histoFactory::fillAttributesData(){
  att[0] = fqData->fq1rnll[0][2]-fqData->fq1rnll[0][1];
  att[1] = fqData->fq1rnll[1][2]-fqData->fq1rnll[1][1];
  return;
}

void histoFactory::fillAttributesMC(){
  att[0] = fqMC->fq1rnll[0][2]-fqMC->fq1rnll[0][1];
  att[1] = fqMC->fq1rnll[1][2]-fqMC->fq1rnll[1][1];
  return;
}



void histoFactory::fillHistos(){
  int nevdata = dataTree->GetEntries();
  int nevmc   = mcTree->GetEntries();
  //fill data histos
  for (int i=0;i<nevdata;i++){
    dataTree->GetEntry(i);
    fillAttributesData();
    for (int iatt=0;iatt<nAttributes;iatt++){
      hManager->fillHistogramData(fqData->nsample,fqData->nbin,iatt,att[iatt],fqData->evtweight);
   }
  }
  //fill MC histos
  for (int j=0;j<nevmc;j++){
    mcTree->GetEntry(j);
    fillAttributesMC();
    for (int jatt=0;jatt<nAttributes;jatt++){
      hManager->fillHistogram(fqMC->nsample,fqMC->nbin,fqMC->ncomponent,jatt,att[jatt],fqMC->evtweight);
    }
  } 
  return;
}

void histoFactory::setDataTree(TChain* ch){
  dataTree=(TTree*)ch;
  fqData = new fQreader(dataTree);
  return;
}

void histoFactory::setDataTree(TTree* tr){
  dataTree=tr;
  fqData = new fQreader(dataTree);
  return;
}

void histoFactory::setMCTree(TTree* tr){
  mcTree=tr;
  fqMC = new fQreader(mcTree);
  return;
}

void histoFactory::setMCTree(TChain* ch){
  mcTree=(TTree*)ch;
  fqMC = new fQreader(mcTree);
  return;
}

//writes all histograms to output file
void histoFactory::saveToFile(){
  fout->Write();
  return;
}

