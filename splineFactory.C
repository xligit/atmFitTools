#include "hSplines.C"
//#include "shared.h"
//#include "fQreader.C"
#include "histoManager.C"
#include "sharedPars.C"
#include "atmFitPars.C"

#define NSYSTMAX 10
#define NHBINSMAX 300
#define NPTSMAX   5

using namespace std;

//class for creating splines
class splineFactory{
  public:
  //constructor
  splineFactory(int nsamp, int nbin, int ncomp, int natt, int nsyst, const char* name);
  splineFactory(const char* parfile);//< initialize using parameters in par file
  splineFactory(){;}
  
  //internal vars
  TString parFileName; //< name of parameter file
  TString nameTag; //< set in constructor. this is the prefix for the output file
  TTree* mcTree; 
  fQreader* mcEvt;
  histoManager* hManager; //manages all default histograms
  TH1D* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; //array for modified histograms for spline creation
  void resetModHistos();
  TFile *fout; //output file
  TString foutName;
  TH1D* htmp; //temporary histogram pointer
  int nSamp;
  int nBin;
  int nComp;
  int nAtt;
  int nSyst;
  TString sysParType; //< code denoting the type of parameterization used, see setupSysPar
  double sysPar[NSYSTMAX]; //systematic parameter values
  double sysUnc[NSYSTMAX];  //systematic parameter uncertainties
  atmFitPars* fitPars;
  double attribute[NATTMAX];
  double eventWeight;
  sharedPars* runpars; //< runtime parameters

  //for output tree
  TTree* splineTree;
  int nbin;
  int ncomponent;
  int nattribute;
  int nsample;
  int nsystpar;
  int npoints;
  int nhistobins;
  double systParValues[NPTSMAX];
  double binWeight[NPTSMAX][NHBINSMAX];
  //methods
  //this needs to be modified for each systematic paramater to add
  double getEvtWeight(int ipar); //returns event weight after applying syst. par. 
  double getEvtWeight(fQreader* mcEvt,int ipar,double value); //
  void setOutputFileName(const char* name){foutName=name;}
  TString getOutputFileName(){return foutName;}
  //
  void fillHistograms(int ipt,int isyst); //fills all histograms given weight
  void  makeManagerFromFile(const char* fname); //reads in histograms from histoFactory
  void fillLeaves(int nsamp,int nbin,int ncomp,int natt,int isyst); //fills leaves of output tree
  void setMCTree(TTree* tr);
  //build the splines
  void buildTheSplines();

  //debugging
  void debugtest();

  //do everything
  void runSplineFactory();

//  private:
  void fillAttributes();
  void setupHistos();
  void setupSystPars(); //sets up systematic parameters
  void incrementSystPars(double nsig);
};

////////////////////////////////////////////////
//run spline factory using parameters from file
void splineFactory::runSplineFactory(){
  
  //create histogram manager from file with prefilled histograms
  makeManagerFromFile(runpars->hFactoryOutput.Data());

  //initialize histogram arrays using the hManager as a template
  setupHistos();

  //Initialize systematic parameter values
 // setupSystPars();
  
  //get pointer to MC tree
  TChain chmc("h1");
  chmc.Add(runpars->hFactoryMCFiles.Data());
  setMCTree((TTree*)&chmc);
  
  //make them splines
  buildTheSplines();

  ////////////////////
  return;
}


void splineFactory::resetModHistos(){
  for (int ibin=0;ibin<nBin;ibin++){
    for (int isamp=0;isamp<nSamp;isamp++){
      for (int icomp=0;icomp<nComp;icomp++){
        for (int iatt=0;iatt<nAtt;iatt++){
          for (int ipt=0;ipt<NPTSMAX;ipt++){
            hMC[isamp][ibin][icomp][iatt][ipt]->Reset();
          }
        }
      }
    }
  }
  return;
}

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
   //fills branches in spline information tree
   nsample = isamp;
   nbin = ibin;
   ncomponent=icomp;
   nattribute=iatt;
   nhistobins = hMC[isamp][ibin][icomp][iatt][0]->GetNbinsX();
   npoints = NPTSMAX;
   nsystpar = isyst;
   for (int ipt=0;ipt<NPTSMAX;ipt++){
     for (int jhistobin=0;jhistobin<=nhistobins;jhistobin++){
       if (hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin)==0){
         binWeight[ipt][jhistobin] = 1.;
       }
       else{
         binWeight[ipt][jhistobin] =
           ((double)hMC[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
           (double)hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin);
       }
     }
   }
   //fill array of all systematic parameter points used to fill histograms
   double sigvals[5] = {-5.,-2.,0.,2.,5.};
   for (int jpt=0;jpt<NPTSMAX;jpt++){
     incrementSystPars(sigvals[jpt]);
  //   cout<<"par: "<<isyst<<sysPar[isyst]<<endl;
     systParValues[jpt]=sysPar[isyst];
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build splines using event-by-event reweighting for various error 
void splineFactory::buildTheSplines(){

  // file setup
  TString fname = nameTag.Data();
  fname.Append("_splineOut.root");
  fout = new TFile(foutName.Data(),"recreate");

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
  splineTree->Branch("systParValues",systParValues,Form("systParValues[%d]/D",NPTSMAX));
  splineTree->Branch("binWeight",binWeight,Form("binWeight[%d][%d]/D",NPTSMAX,NHBINSMAX));
  

  //setup systematic deviations (in sigma)
  double sigvals[5] = {-4.,-2.,0.,2.,4.};

  cout<<"creating spines"<<endl; 

  //loop over each systematic parameter specified in the parameter file
  for (int isyst=0;isyst<fitPars->nSysPars;isyst++){
    resetModHistos(); //< sets all bin contents to zero
    //loop over the MC events
    for (int iev=0;iev<mcTree->GetEntries();iev++){
      mcTree->GetEvent(iev);  //read event
      fillAttributes();  //fill all attributes of this event
      //fill histogram for each point
      for (int ipt=0;ipt<NPTSMAX;ipt++){
        incrementSystPars(sigvals[ipt]);
        getEvtWeight(isyst); 
        //fill histogram fo each attribute
        for (int iatt=0;iatt<nAtt;iatt++){
          hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
        } //end attribute loop
      }  //end point loop
    }//end event loop
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

void splineFactory::incrementSystPars(double nsig){

  //reset initial parameters
  setupSystPars();

  //change systematic parameters
  for (int isyst=0;isyst<nSyst;isyst++){
    sysPar[isyst] += sysUnc[isyst]*nsig;
   // if (sysPar[isyst]<0.) sysPar[isyst]=0.;
  }
  
  return;
}

void splineFactory::setupSystPars(){

  if (!sysParType.CompareTo("tn186")){
    nSyst=0;
    //CCQE xsec norm bin 1//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 1.0;
    nSyst++;
    //CCQE xsec norm  bin 2//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.25;
    nSyst++;
    //CCQE xsec norm bin 3//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.1;
    nSyst++;
    //CCQE xsec norm bin 4//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.05;
    nSyst++;
    //SubGeV flux norm//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.25;
    nSyst++;
    //MultiGeV flux norm//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.15;
    nSyst++;
    //CCnQE xsec norm//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.2;
    nSyst++;
    //NC xsec norm
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.2;
    nSyst++;
    //mu/e xsec ratio
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.05;
    nSyst++;
  }

    return;
}


///////////////////////////////////////////////////////////////////////
//Initialize the histogram arrays needed to builde the splines
void splineFactory::setupHistos(){
  cout<<"called setupHistos()"<<endl;
  TH1D* htemplate;
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
            hMC[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
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
  
  //reset bin contents
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jsamp=0;jsamp<nSamp;jsamp++){
      for (int jatt=0;jatt<nAtt;jatt++){
        for (int jcomp=0;jcomp<nComp;jcomp++){
          hMC[jsamp][jbin][jcomp][jatt][ipt]->Reset();
        }
      }
    }
  }

  //fill new bin contents
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


/////////////////////////////////////////////////////////////////////////////////////////////
//Re-weight the event given the value "value" of additional systematic parameter "ipar"
//Character keys determine the type of systematic error parameters to use
double splineFactory::getEvtWeight(fQreader* mcEvt,int ipar,double value){
  
  double ww = mcEvt->evtweight;
  int absmode = TMath::Abs(mcEvt->mode);
  double Enu     = mcEvt->pmomv[0];
  int  nutype  = TMath::Abs(mcEvt->ipnu[0]);  
 
  //simple TN186 parameterization
  if (!sysParType.CompareTo("tn186")){ 
    //CCQE norm bin1 
    if (ipar==0){
      if ((absmode==1)&&(Enu<200.)) ww*=value;
    }
    //CCQE norm bin2 
    if (ipar==1){
      if ((absmode==1)&&(Enu>200.)&&(Enu<400.)) ww*=value;
    }
    //CCQE norm bin3 
    if (ipar==2){
      if ((absmode==1)&&(Enu>400.)&&(Enu<800.)) ww*=value;
    }
    //CCQE norm bin4 
    if (ipar==3){
      if ((absmode==1)&&(Enu>800.)) ww*=value;
    }
    //SubGevFlux
    if (ipar==4){
      if (Enu<1000.) ww*=value;
    }
    //MultiGeVFlux
    if (ipar==5){
      if (Enu>1000.) ww*=value;
    }
    //CCnQE
    if (ipar==6){
      if ((absmode>1)&&(absmode<30)) ww*=value;
    }
    //NC
    if (ipar==7){
      if (absmode>=30) ww*=value;
    }
    //mu2e ratio
    if (ipar==8){
      if (nutype==14) ww*=value;
    }
  }
  if (ww<0.) ww = 0.;
  eventWeight = ww;
  return ww;

}

double splineFactory::getEvtWeight(int ipar){
//  double ww = 1.;
  double ww = mcEvt->evtweight;
  int absmode = TMath::Abs(mcEvt->mode);
  double Enu     = mcEvt->pmomv[0];
  int  nutype  = TMath::Abs(mcEvt->ipnu[0]);
//  if (ipar==0){
//   if (mcEvt->ncomponent==0) ww*=sysPar[0];
//  } 
//  if (ipar==1){
//   if (mcEvt->ncomponent==1) ww*=sysPar[1];
//  }
//  if (ipar==2){
//   if (mcEvt->ncomponent==2) ww*=sysPar[2];
//  }
   
  //CCQE norm bin1 
  if (ipar==0){
    if ((absmode==1)&&(Enu<200.)) ww*=sysPar[0];
  }
  //CCQE norm bin2 
  if (ipar==1){
    if ((absmode==1)&&(Enu>200.)&&(Enu<400.)) ww*=sysPar[1];
  }
  //CCQE norm bin3 
  if (ipar==2){
    if ((absmode==1)&&(Enu>400.)&&(Enu<800.)) ww*=sysPar[2];
  }
  //CCQE norm bin4 
  if (ipar==3){
    if ((absmode==1)&&(Enu>800.)) ww*=sysPar[3];
  }
  //SubGevFlux
  if (ipar==4){
    if (Enu<1000.) ww*=sysPar[4];
  }
  //MultiGeVFlux
  if (ipar==5){
    if (Enu>1000.) ww*=sysPar[5];
  }
  //CCnQE
  if (ipar==6){
    if ((absmode>1)&&(absmode<30)) ww*=sysPar[6];
  }
  //NC
  if (ipar==7){
    if (absmode>=30) ww*=sysPar[7];
  }
  //mu2e ratio
  if (ipar==8){
    if (nutype==14) ww*=sysPar[8];
  }

  if (ww<0.) ww = 0.;
  eventWeight = ww;
  return ww;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Construct from parameter file
splineFactory::splineFactory(const char*  parfile){

  //set file name and read in pars
  parFileName = parfile;
  cout<<"splineFactory: Parameter File: "<<parfile<<endl;
  runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();

  //setup systematic parameters object
  fitPars = new atmFitPars(parfile);

  //fix internal variables
  nSamp=runpars->nSamples;
  nBin=runpars->nFVBins;
  nComp=runpars->nComponents;
  nAtt=runpars->nAttributes;
  nSyst=runpars->nSysPars;
  nameTag=runpars->globalRootName.Data();
  foutName = runpars->splineFactoryOutput.Data(); //< name of output file with spline parameters
  sysParType = runpars->sysParType; 

  return;

}

splineFactory::splineFactory(int isamp, int ibin, int icomp, int iatt, int isyst, const char* name){
  nameTag = name;
  nSamp = isamp;
  nBin  = ibin;
  nComp = icomp;
  nAtt  = iatt;
  nSyst = isyst;
  foutName = "splineFactoryOutput.root";
}

