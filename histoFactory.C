#include "histoFactory.h"


/////////////////////////////////////////////////
//construct from parameter file
histoFactory::histoFactory(const char* parfile){

  //read in parameters  
  runpars = new sharedPars(parfile);
  runpars->readParsFromFile();
  parFileName = parfile;

  //setup factory
  nameTag = runpars->globalRootName;
 // nameTag.Append("_histograms.root");
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

  cout<<"histoFactory: Initializing histograms..."<<endl;  
  //initialize histograms
  init();

  //loop over events and fill histograms
  cout<<"histoFactory: Filling all histograms"<<endl;
  fillHistos();

  //normalize all of the histograms
  normalizeHistos();

  cout<<"histoFactory: Saving histograms"<<endl;
  //save all filled histograms
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

//void histoFactory::addAttribute(int iatt){
  //adds an attribute to the attribute list
  //the attribute type is a code for the fiTQun output being used
  //this code is used to determine histogram binning and names

  //list of codes:
  //0 - e/mu likelihood ratio for subev 1
  //1 - e/mu likelihood ratio for subev 2
  //

 // attType[nAttributes]=iatt; 
 // nAttributes++;
 // return;
//}

/////////////////////////////////////////////////////////////////////////////
//builds all new histograms for data and MC. Call after all attributes are specified
void histoFactory::init(){

  //////////////////////////////////////
  //setup name for output file 
//  if (!outputFileName.CompareTo("")){
//    outputFileName = "histoFactoryOutput.root";
//  }
//  outputFileName = runpars->hFactoryOuputDir.Data();
  outputFileName = runpars->hFactoryOutput.Data();
 // outputFileName.Append( nameTag.Data() );
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
         //cout<<"Making histogram: "<<hname.Data()<<endl;
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
         //  cout<<"Making histogram: "<<hname.Data()<<endl;
           hManager->setHistogram(isamp,ibin,icomp,iatt,0,getHistogram(iatt,hname.Data()));
        }
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////
//use this function to build histograms
TH1D* histoFactory::getHistogramData(int iatt, const char* thename){
 
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

//  int nBinsNllEMu = 50;
//  if (iatt==0){
//     hnew = new TH1D(thename,thename,nBinsNllEMu,-3000,6000);
//  }
//  if (iatt==1){
//     hnew = new TH1D(thename,thename,nBinsNllEMu,-3000,6000);
//  }


  return hnew;
}



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

//  int nBinsNllEMu = 50;
//  if (iatt==0){
//     hnew = new TH1D(thename,thename,nBinsNllEMu,-3000,6000);
//  }
//  if (iatt==1){
//     hnew = new TH1D(thename,thename,nBinsNllEMu,-3000,6000);
//  }


  return hnew;
 // int nBinsNllEMu = 50;
 // int nBinsNllEMuData = 50;
 // if (iatt==0){
 //    hnew = new TH1D(thename,thename,nBinsNllEMu,-3000,6000);
 // }
 // if (iatt==1){
 //    hnew = new TH1D(thename,thename,nBinsNllEMu,-3000,6000);
 // }
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

void histoFactory::normalizeHistos(double scale){
  //scale all MC histograms by some scaling factor
  if (scale < 0.){
    scale = (double)nDataEvents/(double)nMCEvents;
  }
  hnorm = new TH1D("hnorm","hnorm",1,0,1);
  scale = (double)nDataEvents/(double)nMCEvents;
  hnorm->SetBinContent(1,scale);
  for (int ibin=0;ibin<nBins;ibin++){
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        for (int icomp=0;icomp<nComponents;icomp++){
 //         hManager->getHistogram(isamp,ibin,icomp,iatt)->Scale(scale);
        }
      }
    }
  }
  return;
}

void histoFactory::fillHistos(){
  int nevdata = dataTree->GetEntries();
  int nevmc   = mcTree->GetEntries();
  cout<<"histoFactory: Number of Data entries: "<<nevdata<<endl;
  cout<<"histoFactory: Number of MC entries: "<<nevmc<<endl;
  //fill data histos
  for (int i=0;i<nevdata;i++){
    dataTree->GetEntry(i);
    //fillAttributesData();
    for (int iatt=0;iatt<nAttributes;iatt++){
      hManager->fillHistogramData(fqData->nsample,fqData->nbin,iatt,
                                  fqData->attribute[iatt],fqData->evtweight);
   }
  }
  //fill MC histos
  for (int j=0;j<nevmc;j++){
    mcTree->GetEntry(j);
    fillAttributesMC();
    for (int jatt=0;jatt<nAttributes;jatt++){
      hManager->fillHistogram(fqMC->nsample,fqMC->nbin,fqMC->ncomponent,
                              jatt,fqMC->attribute[jatt],fqMC->evtweight);
    //  hManager->fillHistogram(fqMC->nsample,fqMC->nbin,fqMC->ncomponent,jatt,att[jatt],fqMC->evtweight);
    }
  } 
  return;
}

void histoFactory::setDataTree(TChain* ch){
  dataTree=(TTree*)ch;
  fqData = new fQreader(dataTree);
  nDataEvents = dataTree->GetEntries();
  return;
}

void histoFactory::setDataTree(TTree* tr){
  dataTree=tr;
  fqData = new fQreader(dataTree);
  nDataEvents = dataTree->GetEntries();
  return;
}

void histoFactory::setMCTree(TTree* tr){
  mcTree=tr;
  fqMC = new fQreader(mcTree);
  nMCEvents = mcTree->GetEntries();
  return;
}

void histoFactory::setMCTree(TChain* ch){
  mcTree=(TTree*)ch;
  fqMC = new fQreader(mcTree);
  nMCEvents = mcTree->GetEntries();
  return;
}

//writes all histograms to output file
void histoFactory::saveToFile(){
  cout<<"histoFactory: Saving file: "<<fout->GetPath()<<endl;
  fout->Write();
  return;
}

