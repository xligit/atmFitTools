#ifndef SPLINEFACTORY_C
#define SPLINEFACTORY_C

#include "splineFactory.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//run spline factory using parameters from file
void splineFactory::runSplineFactory(){
  

  
  //make them splines
  buildTheSplines();

  ////////////////////
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::setMCTree(TTree* tr){
  mcTree = tr;
  mcEvt = new fqProcessedEvent(mcTree);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::debugtest(){
  //create new factory
  splineFactory* sfactory = new splineFactory(3,3,7,1,1,"debugtest");
  //set histogram manager and create template histograms;
  sfactory->makeManagerFromFile("factoryOut_factorytest.root");
  

  return;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fills branches of output tree with information required to build a spline 
void splineFactory::fillBranches(int isamp,int ibin,int icomp,int iatt,int isyst){

   // these numbers identify which histogram bin the spline is for
   nsample = isamp;
   nbin = ibin;
   ncomponent = icomp;
   nattribute = iatt;
   nsystpar = isyst;

   // loop over the evaluation points for this bin to create array of bin weights
   nhistobins = hMC[isamp][ibin][icomp][iatt][0]->GetNbinsX();
   npoints = NPTSMAX;
   for (int ipt=0;ipt<NPTSMAX;ipt++){
     for (int jhistobin=1;jhistobin<=nhistobins;jhistobin++){
       if (hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin)==0){
         binWeight[ipt][jhistobin] = 0.;
       }
       else{
         binWeight[ipt][jhistobin] =
           ((double)hMC[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
           (double)hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin);
       }
     }
   }

   // fill array of all systematic parameter points used to fill histograms
   for (int jpt=0;jpt<NPTSMAX;jpt++){
     incrementSystPars(isyst,sigmaValues[jpt]);
     systParValues[jpt]=fitPars->sysPar[isyst];
   }

   return;
}

void splineFactory::setupOutputFile(){

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
 // splineTree->Branch("h2DWeights","TH2D",&h2DWeights,320000,0);;

 return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build splines using event-by-event reweighting for various error 
void splineFactory::buildTheSplines(){



  // setup output file and initialize output TTree
  setupOutputFile();


  //loop over each systematic parameter specified in the parameter file
  cout<<"creating spines"<<endl; 
  int N = fitPars->nSysPars;
  for (int isyst=0;isyst<N;isyst++){

    // loop over MC events, fill histograms, and save bin weights to TTree 
    buildSplineForPar(isyst);

  } 
 
  //write output tree with spline parameters
  splineTree->Write();

  //close output file
  fout->Close();

  ///////
  return;
}


////////////////////////////////////////////////////////////////////////////////
// Build splines while modifying an individual systematic parameter
void splineFactory::buildSplineForPar(int isyspar){

    // set all bin contents to zer
    resetModHistos(); 

    // talk about it
    cout<<"splineFactory: Filling histograms while modifying systematic parameter: "<<isyspar<<endl;

    //loop over the MC events and re-fill all histograms
    int Nevts = mcTree->GetEntries();
    for (int iev=0;iev<Nevts;iev++){

      // talk about it
      if ((iev%5000)==0) cout<<" splineFactory: Reading MC event "<<iev<<endl;
      
      // read in MC info
      mcTree->GetEvent(iev);  //< read event

      //loop over the values of this parameter
      for (int ipt=0;ipt<NPTSMAX;ipt++){

        // change the value of this parameter by some # of sigmas
        incrementSystPars(isyspar,sigmaValues[ipt]); 

        // get the total event weight including systematic
        getEvtWeight(mcEvt,isyspar,fitPars->sysPar[isyspar]);

        //loop over fiTQun attributes and fill histogram for each one
        for (int iatt=0;iatt<nAtt;iatt++){
          hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]
            ->Fill(mcEvt->attribute[iatt],eventWeight);
        } 
      }  
    }
 
    //loop over all histograms to create spline information
    for (int kbin=0; kbin<nBin; kbin++){
      for (int kcomp=0; kcomp<nComp; kcomp++){
        for (int ksamp=0; ksamp<nSamp; ksamp++){
          for (int katt=0; katt<nAtt; katt++){
            // fill spline information for this histogram ;
            fillBranches(ksamp,kbin,kcomp,katt,isyspar);
            splineTree->Fill();
          }
        }
      }
    }
    return;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::incrementSystPars(int isyspar, double nsig){

  //reset parameters to defaults
  fitPars->resetDefaults();

  //adjust systematic parameter in question
  double adjustment = nsig*fitPars->sysParUnc[isyspar];
  double value = fitPars->sysPar[isyspar];
  fitPars->setSysParameter(isyspar, (value+adjustment));

  return;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Initialize the histogram arrays needed to build the splines
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
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void splineFactory::makeManagerFromFile(const char* fname){
  hManager = new histoManager(fname,nSamp,nBin,nComp,nAtt);
  return;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Re-weight the event given the value "value" of additional systematic parameter "ipar"
//Wraps indepenantly defined function "getSystWeight"
double splineFactory::getEvtWeight(fqProcessedEvent* mcevent,int ipar,double value){
  
  //first get the original weight of this event:
  double ww = mcevent->evtweight;

  //get the additional weight from systematic parameters
  double systweight = getSystWeight(sysParType.Data(),mcevent,ipar,value);

  //apply additional weight
  ww *= systweight;

  //set internal variable
  eventWeight = ww; 

  //////////
  return ww;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Construct from parameter file
splineFactory::splineFactory(const char*  parfile){

  //set file name and read in pars
  parFileName = parfile;
  cout<<"splineFactory: Parameter File: "<<parfile<<endl;
  runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
;
  //setup systematic parameters object
  fitPars = new atmFitPars(parfile);

  //fix internal variables
  nSamp=runpars->nSamples;
  nBin=runpars->nFVBins;
  nComp=runpars->nComponents;
  nAtt=runpars->nAttributes;
  nSyst=fitPars->nSysPars;
  nameTag=runpars->globalRootName.Data();
  foutName = runpars->splineFactoryOutput.Data(); //< name of output file with spline parameters
  sysParType = runpars->sysParType; 

  //create histogram manager from file with prefilled histograms
  makeManagerFromFile(runpars->hFactoryOutput.Data());

  //initialize histogram arrays using the hManager as a template
  setupHistos();

  //get pointer to MC tree
  TChain *chmc = new TChain("h1");
  chmc->Add(runpars->hFactoryMCFiles.Data());
  setMCTree((TTree*)chmc);

  //define sigma values
  sigmaValues[0] = -4.;
  sigmaValues[1] = -2.;;
  sigmaValues[2] = 0.;
  sigmaValues[3] = 2.;
  sigmaValues[4] = 4.;


  return;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Original constructor
splineFactory::splineFactory(int isamp, int ibin, int icomp, int iatt, int isyst, const char* name){
  nameTag = name;
  nSamp = isamp;
  nBin  = ibin;
  nComp = icomp;
  nAtt  = iatt;
  nSyst = isyst;
  foutName = "splineFactoryOutput.root";
}


#endif
