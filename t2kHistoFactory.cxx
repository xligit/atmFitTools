#include "t2kHistoFactory.h"


/////////////////////////////////////////////////
//construct from parameter file
t2kHistoFactory::t2kHistoFactory(const char* parfile, bool separateneutmode)
  : separateNeutMode(separateneutmode)
{

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
  if (separateNeutMode) {
    nModes = NMODE;
    hManager = new histoManager(nSamples,nBins,nComponents,nameTag.Data(), nModes, separateNeutMode);
  } else {
    nModes = 0;
    hManager = new histoManager(nSamples,nBins,nComponents,nameTag.Data()); 
  }

}

////////////////////////////////////////////////////////////////
//run histogram factory using the parameter file
void t2kHistoFactory::runHistoFactory(){

  //setup trees and chains for data and mc
  cout<<"t2kHistoFactory: Setting up chains for data and MC.."<<endl;
  TChain *chmc = new TChain("h1");
  TChain *chdata = new TChain("h1");
  chmc->Add(runpars->hFactoryMCFiles.Data());
  setMCTree((TChain*)chmc);
  chdata->Add(runpars->hFactoryDataFiles.Data());  
  setDataTree((TChain*)chdata);

  cout<<"t2kHistoFactory: Initializing histograms..."<<endl;  
  //initialize histograms
  init();

  //loop over events and fill histograms
  cout<<"t2kHistoFactory: Filling all histograms"<<endl;
  fillHistos();

  //normalize all of the histograms
  normalizeHistos();

  cout<<"t2kHistoFactory: Saving histograms"<<endl;
  //save all filled histograms
  saveToFile();
  chmc->Delete();
  chdata->Delete();
  ///////////////////
  return;
}


////////////////////////////////////////////////////////////////
//constructor to re-created a histogram factory from a file
t2kHistoFactory::t2kHistoFactory(int nsampl,int nbins,int ncomp,const char* name){
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

t2kHistoFactory::t2kHistoFactory(int nsampl, int nbins, int ncomp, int nmode, bool separateneutmode, const std::string name)
  : nSamples(nsampl)
  , nComponents(ncomp)
  , nBins(nbins)
  , nModes(nmode)
  , separateNeutMode(separateneutmode)
{
  nameTag = name.c_str();
  nameTag.Append("_hFactoryOutPut");
  nAttributes = 0;
  hManager = new histoManager(nsampl, nbins, ncomp, name.c_str(), nmode, separateneutmode);
}
//void t2kHistoFactory::addAttribute(int iatt){
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
void t2kHistoFactory::init(){

  //////////////////////////////////////
  //setup name for output file 
//  if (!outputFileName.CompareTo("")){
//    outputFileName = "t2kHistoFactoryOutput.root";
//  }
//  outputFileName = runpars->hFactoryOuputDir.Data();
  outputFileName = runpars->hFactoryOutput.Data();
 // outputFileName.Append( nameTag.Data() );
  fout = new TFile(outputFileName.Data(),"RECREATE");


  ///////////////////////////////////////////
  //make sure to set sumw2 to proporly calculate errors (?)
  TH1::SetDefaultSumw2(kTRUE);

  if (!separateNeutMode) {
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
  } else {
    // setup data histos
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
	for (int iatt = 0; iatt < nAttributes; ++iatt) {
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
	  for (int imode = 0; imode < nModes; ++imode) {
	    for (int iatt=0;iatt<nAttributes;iatt++){
	      TString hname = "hmc_";
	      hname.Append(Form("samp%d_bin%d_comp%d_mode%d_att%d",isamp,ibin,icomp,imode,iatt));
	      //  cout<<"Making histogram: "<<hname.Data()<<endl;
	      hManager->setHistogram(isamp,ibin,icomp,imode,iatt,0,getHistogram(iatt,hname.Data()));
	      hname.Append("nom");
	      hManager->setNominalHistogram(isamp,ibin,icomp,imode,iatt,getHistogram(iatt,hname.Data()));
	    }
	  }
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////
//use this function to build histograms
TH1D* t2kHistoFactory::getHistogramData(int iatt, const char* thename){
 
  //setup binning
  TString parname = Form("nBinsAtt%d",iatt);
  int nbins = runpars->getParI(parname.Data());
  parname = Form("xMinAtt%d",iatt);
  double xmin = runpars->getParD(parname.Data());
  parname = Form("xMaxAtt%d",iatt);
  double xmax = runpars->getParD(parname.Data());


  //histogram building
  cout<<"t2kHistoFactory: Creating histogram: "<<endl;
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



TH1D* t2kHistoFactory::getHistogram(int iatt, const char* thename){

  //setup binning
  TString parname = Form("nBinsAtt%d",iatt);
  int nbins = runpars->getParI(parname.Data());
  parname = Form("xMinAtt%d",iatt);
  double xmin = runpars->getParD(parname.Data());
  parname = Form("xMaxAtt%d",iatt);
  double xmax = runpars->getParD(parname.Data());


  //histogram building
  cout<<"t2kHistoFactory: Creating histogram: "<<endl;
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
  //return hnew;
}

//calculates attributes from the raw fiTQun output
void t2kHistoFactory::fillAttributesData(){
  att[0] = fqData->fq1rnll[0][2]-fqData->fq1rnll[0][1];
  att[1] = fqData->fq1rnll[1][2]-fqData->fq1rnll[1][1];
  return;
}

void t2kHistoFactory::fillAttributesMC(){
  att[0] = fqMC->fq1rnll[0][2]-fqMC->fq1rnll[0][1];
  att[1] = fqMC->fq1rnll[1][2]-fqMC->fq1rnll[1][1];
  return;
}

void t2kHistoFactory::normalizeHistos(double scale){
  //scale all MC histograms by some scaling factor
  if (scale < 0.){
    scale = (double)nDataEvents/(double)nMCEvents;
  }
  hnorm = new TH1D("hnorm","hnorm",1,0,1);
  scale = (double)nDataEvents/(double)nMCEvents;
  hnorm->SetBinContent(1,scale);
  /*for (int ibin=0;ibin<nBins;ibin++){
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        for (int icomp=0;icomp<nComponents;icomp++){
 //         hManager->getHistogram(isamp,ibin,icomp,iatt)->Scale(scale);
        }
      }
    }
    }*/
  return;
}

void t2kHistoFactory::fillHistos(){
  int nevdata = dataTree->GetEntries();
  int nevmc   = mcTree->GetEntries();
  cout<<"t2kHistoFactory: Number of Data entries: "<<nevdata<<endl;
  cout<<"t2kHistoFactory: Number of MC entries: "<<nevmc<<endl;
  //fill data histos
  for (int i=0;i<nevdata;i++){
    dataTree->GetEntry(i);
    //fillAttributesData();
    if (i%100==0) {
      std::cout<<i<<" "<<fqData->nsample<<" "<<fqData->nbin<<" "<<fqData->attribute[0]<<" "<<fqData->attribute[1]<<std::endl;
    }
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
      if (!separateNeutMode) {
	hManager->fillHistogram(fqMC->nsample,fqMC->nbin,fqMC->ncomponent,
				jatt,fqMC->attribute[jatt],fqMC->rfgweight);
      } else {
	hManager->fillHistogram(fqMC->nsample, fqMC->nbin, fqMC->ncomponent, fqMC->nmode, jatt, fqMC->attribute[jatt], fqMC->rfgweight);
	hManager->fillNominalHistogram(fqMC->nsample, fqMC->nbin, fqMC->ncomponent, fqMC->nmode, jatt, fqMC->attribute[jatt], fqMC->evtweight);
      }
    //  hManager->fillHistogram(fqMC->nsample,fqMC->nbin,fqMC->ncomponent,jatt,att[jatt],fqMC->evtweight);
    }
  } 
  return;
}

void t2kHistoFactory::setDataTree(TChain* ch){
  dataTree=ch;
  fqData = new t2kfqReader(dataTree);
  nDataEvents = dataTree->GetEntries();
  dataTree->GetEntry(0);
  std::cout<<"got data entry 0"<<std::endl;
  return;
}

void t2kHistoFactory::setDataTree(TTree* tr){
  dataTree=(TChain*)tr;
  fqData = new t2kfqReader(dataTree);
  nDataEvents = dataTree->GetEntries();
  return;
}

void t2kHistoFactory::setMCTree(TTree* tr){
  mcTree=(TChain*)tr;
  fqMC = new t2kfqReader(mcTree);
  fqMC->FillMap();
  nMCEvents = mcTree->GetEntries();
  return;
}

void t2kHistoFactory::setMCTree(TChain* ch){
  mcTree=ch;
  fqMC = new t2kfqReader(mcTree);
  fqMC->FillMap();
  nMCEvents = mcTree->GetEntries();
  return;
}

//writes all histograms to output file
void t2kHistoFactory::saveToFile(){
  cout<<"t2kHistoFactory: Saving file: "<<fout->GetPath()<<endl;
  fout->Write();
  return;
}

