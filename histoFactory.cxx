#ifndef HISTOFACTORY_C
#define HISTOFACTORY_C


#include "histoFactory.h"
;

/////////////////////////////////////////////////
//construct from parameter file
histoFactory::histoFactory(const char* parfile){

  //read in parameters  
  runpars = new sharedPars(parfile);
  runpars->readParsFromFile();
  parFileName = parfile;

  //setup factory
  nameTag = runpars->globalRootName;
  nSamples = runpars->nSamples;
  nComponents = runpars->nComponents;
  nAttributes = runpars->nAttributes;
  nBins = runpars->nFVBins;

  //histogram manager setup
  hManager = new histoManager(nSamples,nBins,nComponents,nameTag.Data()); 

}

////////////////////////////////////////////////////////////////
//run histogram factory using the parameter file
void histoFactory::runHistoFactory(){

  //setup trees and chains for data and mc
  cout<<"histoFactory: Setting up chains for data and MC.."<<endl;
  TChain chmc("h1");
  TChain chdata("h1");
  chmc.Add(runpars->hFactoryMCFiles.Data());
  setMCTree((TTree*)&chmc);
  chdata.Add(runpars->hFactoryDataFiles.Data());  
  setDataTree((TTree*)&chdata);
  cout<<"histoFactory: Number of data events: "<<dataTree->GetEntries()<<endl;
  cout<<"histoFactory: Number of MC events: "<<mcTree->GetEntries()<<endl;

  cout<<"histoFactory: Initializing histograms..."<<endl;  

  //Initialize  histograms
  init();

  //loop over events and fill histograms
  cout<<"histoFactory: Filling all histograms"<<endl;
  fillHistos();

  //normalize all of the histograms
  normalizeHistos();

  //save all filled histograms
  cout<<"histoFactory: Saving histograms"<<endl;
  saveToFile();

  ///////////////////
  return;
}


////////////////////////////////////////////////////////////////
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



/////////////////////////////////////////////////////////////////////////////
//builds all new histograms for data and MC. 
void histoFactory::init(){

  //////////////////////////////////////
  //setup name for output file 
  outputFileName = runpars->hFactoryOutput.Data();
  fout = new TFile(outputFileName.Data(),"RECREATE");

  ///////////////////////////////////////////
  //make sure to set sumw2 to proporly calculate errors (?)
  TH1D* hsetsum = new TH1D();
  hsetsum->SetDefaultSumw2(kTRUE);

  ////////////////////////////////////////
  //setup data histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int iatt=0;iatt<nAttributes;iatt++){
         TString hname = "hdata_";
         hname.Append(Form("samp%d_bin%d_att%d",isamp,ibin,iatt));
         hManager->setHistogram(isamp,ibin,0,iatt,1,getHistogramData(iatt,hname.Data()));
      }
    }
  }

  //setup mc histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int icomp=0;icomp<nComponents;icomp++){
        for (int iatt=0;iatt<nAttributes;iatt++){
           TString hname = "hmc_";
           hname.Append(Form("samp%d_bin%d_comp%d_att%d",isamp,ibin,icomp,iatt));
           hManager->setHistogram(isamp,ibin,icomp,iatt,0,getHistogram(iatt,hname.Data()));
        }
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////
//Use this function to build data histograms
TH1D* histoFactory::getHistogramData(int iatt, const char* thename){
 
  //setup binning from runtime parameters
  TString parname = Form("nBinsDataAtt%d",iatt);
  int nbins = runpars->getParI(parname.Data());
  parname = Form("xMinAtt%d",iatt);
  double xmin = runpars->getParD(parname.Data());
  parname = Form("xMaxAtt%d",iatt);
  double xmax = runpars->getParD(parname.Data());


  //histogram building
  cout<<"histoFactory: Creating histogram: "<<endl;
  cout<<"    name:  "<<thename<<endl;
  cout<<"    nbins: "<<nbins<<endl;
  cout<<"    xmin:  "<<xmin<<endl;
  cout<<"    xmax:  "<<xmax<<endl;
  TH1D* hnew = new TH1D(thename,thename,nbins,xmin,xmax);

  //
  return hnew;
}


/////////////////////////////////////////////////////////////////////////////////////
//Use this function to build MC histograms
TH1D* histoFactory::getHistogram(int iatt, const char* thename){

  //setup binning
  TString parname = Form("nBinsAtt%d",iatt);
  int nbins = runpars->getParI(parname.Data());
  parname = Form("xMinAtt%d",iatt);
  double xmin = runpars->getParD(parname.Data());
  parname = Form("xMaxAtt%d",iatt);
  double xmax = runpars->getParD(parname.Data());


  //histogram building
  cout<<"histoFactory: Creating histogram: "<<endl;
  cout<<"    name:  "<<thename<<endl;
  cout<<"    nbins: "<<nbins<<endl;
  cout<<"    xmin:  "<<xmin<<endl;
  cout<<"    xmax:  "<<xmax<<endl;
  TH1D* hnew = new TH1D(thename,thename,nbins,xmin,xmax);

  return hnew;
}

/////////////////////////////////////////////////////
// Find the naive normaliztion
void histoFactory::normalizeHistos(double scale){

  //scale all MC histograms by some scaling factor
  if (scale < 0.){
    cout<<"histoFactory: Finding MC normalization"<<endl;
    double mcsumweights = 0.;
    double datasum = (double)nDataEvents;
    for (int ievt=0;ievt<nMCEvents;ievt++){
      mcTree->GetEntry(ievt);
      mcsumweights += fqMC->evtweight; 
    }
    double datsumweights = 0;
    for (int iev=0; iev<nDataEvents; iev++){
    //  cout<<"evt: "<<iev<<endl;
      dataTree->GetEntry(iev);
    //  cout<<"wgt: "<<fqData->evtweight;
      datsumweights += fqData->evtweight;
    }
    scale = datsumweights/mcsumweights;
    cout<<"   ..norm = "<<datsumweights<<"/"<<mcsumweights<<" = "<<scale<<endl;
  }
  
  hnorm = new TH1D("hnorm","hnorm",1,0,1);
  hnorm->SetBinContent(1,scale);
  return;
}


///////////////////////////////////////////////////////////////////////
// This will fill all of the histrograms built by init()
void histoFactory::fillHistos(){

  // get total numbers of data and MC events
  int nevdata = dataTree->GetEntries();
  int nevmc   = mcTree->GetEntries();
  cout<<"histoFactory: Number of Data entries: "<<nevdata<<endl;
  cout<<"histoFactory: Number of MC entries: "<<nevmc<<endl;

  //fill data histos
  for (int i=0;i<nevdata;i++){
    dataTree->GetEntry(i);
    for (int iatt=0;iatt<nAttributes;iatt++){
      hManager->fillHistogramData(fqData->nsample,fqData->nbin,iatt,
                                  fqData->attribute[iatt],fqData->evtweight);
   }
  }

  //fill MC histos
  for (int j=0;j<nevmc;j++){
    mcTree->GetEntry(j);
    //fillAttributesMC();
    for (int jatt=0;jatt<nAttributes;jatt++){
      hManager->fillHistogram(fqMC->nsample,fqMC->nbin,fqMC->ncomponent,
                              jatt,fqMC->attribute[jatt],fqMC->evtweight);
    }
  } 

  //
  return;
}

void histoFactory::setDataTree(TChain* ch){
  dataTree=(TTree*)ch;
  fqData = new fqProcessedEvent(dataTree);
  nDataEvents = dataTree->GetEntries();
  return;
}

void histoFactory::setDataTree(TTree* tr){
  dataTree=tr;
  fqData = new fqProcessedEvent(dataTree);
  nDataEvents = dataTree->GetEntries();
  return;
}

void histoFactory::setMCTree(TTree* tr){
  mcTree=tr;
  fqMC = new fqProcessedEvent(mcTree);
  nMCEvents = mcTree->GetEntries();
  return;
}

void histoFactory::setMCTree(TChain* ch){
  mcTree=(TTree*)ch;
  fqMC = new fqProcessedEvent(mcTree);
  nMCEvents = mcTree->GetEntries();
  return;
}

//writes all histograms to output file
void histoFactory::saveToFile(){
  cout<<"histoFactory: Saving file: "<<fout->GetPath()<<endl;
  fout->Write();
  return;
}

#endif
