//#include "TTree.h"
//#include "TFile.h"
//#include "TString.h"
//#include "TH1D.h"

//#include "shared.h"

#include "atmFitPars.C"
#include "histoManager.C"
//#include "fQreader.C"
#include "splineFactory.C"
#include "mcmcReader.C"
#include "THStack.h"
#include "sharedPars.C"


#define NCLASSMAX 10

//This class aims to take the MCMC output and generate uncertainties in number of events 
class postfitCalculator{
  public:
  ///////////////////////////////////////////////////////////////////////////////
  //constructor
  //initialize using the file name that the contains the mcmc tree to be analyzed
  //the optional parameter is the name of the shared paremeter file for creating the
  //atmFitPars object.
  postfitCalculator(const char* mcmcfilename,const char* parfile);
  postfitCalculator(const char* parfile); //< construct from parameter file
  void init(); //< called by constructor to initialize histograms, trees, etc

  /////////////////////////////////////////////////////////////
  //internal variables
  TTree* mcmcpath; //< points to tree containing the mcmc steps
  TTree* mctree;  //< points to the tree containg the mc events
  TTree* datatree; //< for data events
  fQreader* datareader; //< points t
  mcmcReader* path; //< points to a step in mcmcpath
  fQreader* mcreader; //< points to an event in mctree
  atmFitPars* fitpars; //< points to the current fit parameters object
  splineFactory* sfact; //< empty spline factory to call getEvtWeight. 
  histoManager* hManager; //< use this for histogram templates and comparing to data
  TFile* outFile; //< outputfile
  int NMCEvents; //< limits # of MC events to use to speed up histo filling
  int NDataEvents; //< limits # of Data events to use to speed up histo filling
  double pars[1000]; 
  double evtweight;
  int evtclass; //< code for event class 
  TString parFileName;  //< name of parameter file
  //for "nsk" tree
  double nevents[10]; //< number of events for each event class
  double neventsdef[10]; //< number of events for each class in default mc
  double neventstot[10]; //< total number of events for each event class
  TTree* nsktree; //< tree to hold # of events for each MCMC point
  int selectiontype; //< code for which selection to use
  int NMCMCPts; //< number of points to sample from
  int MCMCBurnIn; //< number of burn-in points
  sharedPars* runPars;
  int currentMCMCPoint;
  int currentMCEvent;
  double normFactor;
  int nAttributes;
  TCanvas* cc;


  //histogrms of interest
  TH1D* hEnuMu[10][200]; //<energy assuming muon
//  TH1D* hPIDemuSE2[10][200]; //< e/mu PID of second subev
  TH1D* hPIDemu[10][200]; //<e/mu PID
  TH1D* hNSK[10]; //< numbers of events after cuts
  TH1D* hMCMod[3][3][5][20];
  TH1D* hPostFit[3][3][5];
  TH1D* hPassFail[100];
  TH1D* hPassFailData;
  TH1D* hPassFailAvg;
  ////////////////////////////////////////////////////////
  //methods
  // modfies the values of the current MC event using parameters from current MCMC point
  void  initHistos();
  void  modifyCurrentEvent();
  void  setMCTree(TTree* tr);
  void  setDataTree(TTree* tr);
  void  setSelectionType(int itype){selectiontype=itype;}
  double getEvtWeight(); //< returns the weight given the current parameters
  void  modifyAttributes();
  void  unmodifyAttributes();
  void  setParsFromMCMC(int istep); //< sets parameters from mcmmc cloud
  void fillHistos(int iclass,int islot); //< fills all histograms 
  int makeSelection(int iselect,int dataflg =0); //< select events based on various attributes;
  int getEvtClass();
  void drawBreakdown(int ihisto,int islot);
  void makeHistos(int iselection,int islot);
  void findEvtClasses(); //< fills array of total numbers of each event class
  void fillNskTree(int npts,int nburn=1000); //< step through MCMC cloud and fill NSK tree at each point
  void makeNSKTable();
  void showAllHistos(int ihisto,int iclass);
  void MCMCLooperTemplate();
  void cosmicPostFitAnalysis();
  void attributeAnalysis();
  void makeHistoArray(TH1D* harr[], int nhistos, const char* name, int npts, double xmin, double xmax);
  void drawPostFitHisto(int isamp,int ibin, int iatt);
  void drawPostFitHistoAll(int isamp,int ibin, int iatt);
  void calcPostFitHistos(int errtype=0);
  void printPostFitHistos(const char* directory);
//  void runPostFit();
};

void postfitCalculator::printPostFitHistos(const char* directory){

  //loop over histos and pring
  for (int ibin=0;ibin<fitpars->nBins;ibin++){
    for (int isamp=0;isamp<fitpars->nSamples;isamp++){
      for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
         cout<<"break1"<<endl;
         TString plotname = directory;
         cout<<"break2"<<endl;
         plotname.Append("pfit");
         plotname.Append(Form("_postfit_samp%d_bin%d_att%d",isamp,ibin,iatt));
         cout<<"break3"<<endl;
         plotname.Append(".png");
         drawPostFitHisto(isamp,ibin,iatt); //< draw all plot
         cout<<"break4"<<endl;
         cc->Print(plotname.Data());
      }
    }
  }

  return;

}


void postfitCalculator::calcPostFitHistos(int errtype){

  /////////////////////////////////////////
  //initialize hisotgrams
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        TString hname = Form("hMCMCPostfit_samp%d_bin%d_iatt%d",isamp,ibin,iatt);
        hPostFit[isamp][ibin][iatt] = (TH1D*)hMCMod[isamp][ibin][iatt][0]->Clone(hname.Data());
        hPostFit[isamp][ibin][iatt]->SetLineColor(kCyan+1);
        hPostFit[isamp][ibin][iatt]->SetFillColor(kCyan+1);
        hPostFit[isamp][ibin][iatt]->Reset();        
      }
    }
  }

  ////////////////////////////////////////////////////
  //calculate mean 
  double thenorm=(double)NMCMCPts-1.;
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        for (int ipt=1;ipt<NMCMCPts;ipt++){
          //add all other histograms to first one
          hPostFit[isamp][ibin][iatt]->Add(hMCMod[isamp][ibin][iatt][ipt]);
        }
      }
    }
  }
  //now we normalize to get average
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        cout<<"norm: "<<thenorm<<endl;
        hPostFit[isamp][ibin][iatt]->Scale(1./thenorm);
      }
    }
  }

  ////////////////////////////////////////////////////////
  //calculate RMS:
  double rms[3][3][5][100];
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        for (int ipt=1;ipt<NMCMCPts;ipt++){
          for (int hbin=0;hbin<hMCMod[isamp][ibin][iatt][ipt]->GetNbinsX();hbin++){
            //set initial rms to zero
            rms[isamp][ibin][iatt][hbin] = 0.;
          }
        }
      }
    }
  }
  //now we fill rms array
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        for (int ipt=1;ipt<NMCMCPts;ipt++){
          for (int hbin=0;hbin<hMCMod[isamp][ibin][iatt][ipt]->GetNbinsX();hbin++){
            double diff = hMCMod[isamp][ibin][iatt][ipt]->GetBinContent(hbin) - hPostFit[isamp][ibin][iatt]->GetBinContent(hbin);
            rms[isamp][ibin][iatt][hbin] += TMath::Sqrt((diff*diff))/(thenorm-1.);
          }
        }
      }
    }
  }
  //now set error equal to rms
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        for (int hbin=0;hbin<hPostFit[isamp][ibin][iatt]->GetNbinsX();hbin++){
          double staterr = hPostFit[isamp][ibin][iatt]->GetBinError(hbin);
          double systerr = rms[isamp][ibin][iatt][hbin];
          double toterr  = TMath::Sqrt(( staterr*staterr) + (systerr*systerr) );
          if (errtype==1) hPostFit[isamp][ibin][iatt]->SetBinError(hbin,staterr);
          else if (errtype==2) hPostFit[isamp][ibin][iatt]->SetBinError(hbin,systerr);
          else{
            hPostFit[isamp][ibin][iatt]->SetBinError(hbin,toterr);
          }
        }
      }
    }
  }


  return;

}

void postfitCalculator::drawPostFitHisto(int isamp, int ibin, int iatt){
  double norm = (double)datatree->GetEntries()/(double)NMCEvents;
//  hMCMod[isamp][ibin][iatt][0]->SetLineColor(kRed);
//  hMCMod[isamp][ibin][iatt][0]->Scale(norm);
//  for (int ipt=1;ipt<NMCMCPts;ipt++){
 //   hMCMod[isamp][ibin][iatt][ipt]->SetLineColor(kCyan+1);
 //   hMCMod[isamp][ibin][iatt][ipt]->Scale(norm);        
 // }
  hManager->getHistogramData(isamp,ibin,iatt)->SetMarkerStyle(8);   
  hManager->getHistogramData(isamp,ibin,iatt)->SetLineWidth(3);   
  hManager->getHistogramData(isamp,ibin,iatt)->Draw();   
  hPostFit[isamp][ibin][iatt]->Draw("samee2");
  hMCMod[isamp][ibin][iatt][0]->SetMarkerStyle(8);
  hMCMod[isamp][ibin][iatt][0]->SetMarkerColor(kRed);
  hMCMod[isamp][ibin][iatt][0]->SetLineWidth(2);
  hMCMod[isamp][ibin][iatt][0]->Draw("same");
//  hPostFit[isamp][ibin][iatt]->Draw("samee2");
  hManager->getHistogramData(isamp,ibin,iatt)->Draw("same");   

 // for (int ipt=0;ipt<NMCMCPts;ipt++){
 //   hMCMod[isamp][ibin][iatt][ipt]->Scale(1./norm);        
//  }

  return; 
}

void postfitCalculator::drawPostFitHistoAll(int isamp, int ibin, int iatt){
  double norm = (double)datatree->GetEntries()/(double)NMCEvents;
//  hMCMod[isamp][ibin][iatt][0]->SetLineColor(kRed);
//  hMCMod[isamp][ibin][iatt][0]->Scale(norm);
//  for (int ipt=1;ipt<NMCMCPts;ipt++){
 //   hMCMod[isamp][ibin][iatt][ipt]->SetLineColor(kCyan+1);
 //   hMCMod[isamp][ibin][iatt][ipt]->Scale(norm);        
 // }
  hManager->getHistogramData(isamp,ibin,iatt)->SetMarkerStyle(8);   
  hManager->getHistogramData(isamp,ibin,iatt)->SetLineWidth(3);   
  hManager->getHistogramData(isamp,ibin,iatt)->Draw();   

  for (int ipt=1;ipt<NMCMCPts;ipt++){
    hMCMod[isamp][ibin][iatt][ipt]->SetLineWidth(3);        
    hMCMod[isamp][ibin][iatt][ipt]->Draw("samehisto");        
  }
  hMCMod[isamp][ibin][iatt][0]->SetMarkerStyle(8);
  hMCMod[isamp][ibin][iatt][0]->SetMarkerColor(kRed);
  hMCMod[isamp][ibin][iatt][0]->SetLineWidth(2);
  hMCMod[isamp][ibin][iatt][0]->Draw("same");

 // for (int ipt=0;ipt<NMCMCPts;ipt++){
 //   hMCMod[isamp][ibin][iatt][ipt]->Scale(1./norm);        
//  }

  return; 
}

void postfitCalculator::initHistos(){

  //initialize an array of histograms to be filled at each MCMC pt
  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        for (int ipt=0;ipt<NMCMCPts;ipt++){
          TString hname = Form("hMCMCMod_samp%d_bin%d_iatt%d_mcmcpt%d",isamp,ibin,iatt,ipt);
          hMCMod[isamp][ibin][iatt][ipt] = (TH1D*)hManager->getHistogram(isamp,ibin,0,iatt)->Clone(hname.Data());
          hMCMod[isamp][ibin][iatt][ipt]->Reset();
        }
      }
    }
  }

  return;
}

void postfitCalculator::makeHistoArray(TH1D* harr[], int nhistos, const char* name, int npts, double xmin, double xmax){

  
  for (int i=0;i<nhistos;i++){
    TString hname = name;
    hname.Append(Form("_%d",i));
    harr[i] = new TH1D(hname.Data(),hname.Data(),npts,xmin,xmax);
  }

  return;
}

void postfitCalculator::attributeAnalysis(){

  //get number of total steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-MCMCBurnIn;
  if (nsteps<(NMCMCPts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  //loop over points (uniform sampling)
  int dstep = nsteps/NMCMCPts;
  int istep = MCMCBurnIn;
  int samppts[10000];

  //get array of points to sample from
  for (int i=0;i<NMCMCPts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

  //find the number of MC events to use
  if ((NMCEvents<=0)||(NMCEvents>mctree->GetEntries())) NMCEvents = mctree->GetEntries(); 


  //setup histograms to be filled
/*  int    nbinatt0=50;
  double xminatt0=-2000;
  double xmaxatt0=2000;
  int    nbinatt1=50;
  double xminatt1=-1000;
  double xmaxatt1=1000;
  int    nbinatt2=50;
  double xminatt2=-1000;
  double xmaxatt2=1000;
  int    nbinatt3=50;
  double xminatt3=-1000;
  double xmaxatt3=1000;
  const int NHISTO = NMCMCPts;
  TH1D* hatt0[NHISTO];
  TH1D* hatt1[NHISTO];
  TH1D* hatt0data;
  TH1D* hatt1data;
  TString attname0  = runPars->fQAttName0;
  TString attname1  = runPars->fQAttName1;
  makeHistoArray(hatt0,NMCMCPts,attname0.Data(),nbinatt0,xminatt0,xmaxatt0);
  makeHistoArray(hatt1,NMCMCPts,attname1.Data(),nbinatt1,xminatt1,xmaxatt1);
  hatt0data = new TH1D("hatt0data","hatt0data",nbinatt0,xminatt0,xmaxatt0);
  hatt1data = new TH1D("hatt1data","hatt1data",nbinatt1,xminatt1,xmaxatt1); */

  
  //loop for MC defaults
  fitpars->resetDefaults(); //< set parameters to default values (no modifications to MC)
  for (int ievt=0;ievt<NMCEvents;ievt++){
    mctree->GetEntry(ievt);
    for (int iatt=0;iatt<runPars->nAttributes;iatt++){
      hMCMod[mcreader->nsample][mcreader->nbin][iatt][0]->Fill(mcreader->attribute[iatt],mcreader->evtweight);
    }
  }
  //loop over data and fill histos
//  for (int ievt=0;ievt<NDataEvents;ievt++){
//    datatree->GetEvent(ievt);
//    hatt0data->Fill(datareader->attribute[0]);
//    hatt1data->Fill(datareader->attribute[1]);
//  }

  //loop over number of MCMC points
  for (int ipt=1;ipt<NMCMCPts;ipt++){
    currentMCMCPoint = samppts[ipt];
    //loop over MC events
    for (int ievt=0;ievt<NMCEvents;ievt++){
      currentMCEvent = ievt;
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      modifyCurrentEvent();//< all attributes are modified and eventWeight is calculated
      //feel free to do things here, like fill histograms or something  
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        hMCMod[mcreader->nsample][mcreader->nbin][iatt][ipt]->Fill(mcreader->attribute[iatt],evtweight);
      }
    } 
  }


  // calculate normalization of MC using sum of weights
  double mcsumweights = 0.;
  for (int ievt=0;ievt<NMCEvents;ievt++){
    currentMCEvent = ievt;
    mctree->GetEntry(currentMCEvent); //< read in default MC event
    mcsumweights += mcreader->evtweight;
  }

  double norm = (double)datatree->GetEntries()/mcsumweights;

  for (int ibin=0;ibin<runPars->nFVBins;ibin++){
    for (int isamp=0;isamp<runPars->nSamples;isamp++){
      for (int iatt=0;iatt<runPars->nAttributes;iatt++){
        hMCMod[isamp][ibin][iatt][0]->SetLineColor(kRed);
        hMCMod[isamp][ibin][iatt][0]->Scale(norm);
        for (int ipt=1;ipt<NMCMCPts;ipt++){
          hMCMod[isamp][ibin][iatt][ipt]->SetLineColor(kCyan+1);
          hMCMod[isamp][ibin][iatt][ipt]->Scale(norm);        
        }
      }
    }
  }

  //draw it all on same canvas
/*  TCanvas* cc = new TCanvas("cc","cc",800,700);
  hatt0data->Draw();
  hatt0[0]->SetLineColor(kRed);
  hatt0[0]->Draw("same");
  for (int ipt =1;ipt<NMCMCPts;ipt++){
    hatt0[ipt]->SetLineColor(kCyan+1);
    hatt0[ipt]->Draw("same");
  }
  cc->Print("att0.png");
  hatt1data->Draw();
  hatt1[0]->SetLineColor(kRed);
  hatt1[0]->Draw("same");
  for (int ipt =1;ipt<NMCMCPts;ipt++){
    hatt1[ipt]->SetLineColor(kCyan+1);
    hatt1[ipt]->Draw("same");
  }
  cc->Print("att1.png");
*/

  return;

}

///////////////////////////////////////////////////
//runs a post-fit analysis for the stopping cosmics
void  postfitCalculator::cosmicPostFitAnalysis(){

  //////////////////////////////////////////
  //get number of total steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-MCMCBurnIn;
  if (nsteps<(NMCMCPts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  ////////////////////////////////////////
  //setup uniform sampling of MCMC
  int dstep = nsteps/NMCMCPts;
  int istep = MCMCBurnIn;
  int samppts[10000];
  //get array of points to sample from
  for (int i=0;i<NMCMCPts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

  //////////////////////////////////////////////
  //find the number of MC events to use
  if ((NMCEvents<=0)||(NMCEvents>mctree->GetEntries())) NMCEvents = mctree->GetEntries(); 

  ///////////////////////////////////////////////
  //histogram setup
  hPassFailData =  new TH1D("hpassfaildata","hpassfaildata",2,0,1.1);
  hPassFailAvg =  new TH1D("hpassfailavg","hpassfailavg",2,0,1.1);
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    TString hname = "passfail";
    hname.Append("_point%d",ipt);
    hPassFail[ipt] = new TH1D(hname.Data(),hname.Data(),2,0,1.1);
  }

  ///////////////////////////////////////////////////////////
  //loop over number of MCMC points
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    currentMCMCPoint = samppts[ipt];
    //loop over MC events
    for (int ievt=0;ievt<NMCEvents;ievt++){
      currentMCEvent = ievt;
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      modifyCurrentEvent();//< all attributes are modified and eventWeight is calculated
      if (makeSelection(1)){
        hPassFail[ipt]->Fill(1);
      }
      else{
        hPassFail[ipt]->Fill(0);
      }
    } 
  }

  ////////////////////////////////////////////////////////////////
  //loop over data
  for (int ievt=0;ievt<NDataEvents;ievt++){
    datatree->GetEntry(ievt);
    if (makeSelection(1,1)){
       hPassFailData->Fill(1);
    }
    else{
      hPassFailData->Fill(0);
    }
    //feel free to do things here, like fill histograms or something   
  } 

  //////////////////////////////////////////////////////
  //normalize MC
  double norm = (double)NDataEvents/(double)NMCEvents; 
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    hPassFail[ipt]->Scale(norm);
    hPassFail[ipt]->SetLineColor(kBlue);
    hPassFail[ipt]->SetLineWidth(2);
  }

  /////////////////////////////////////////////////////////
  //calculate average and variance over MCMC points
  double mean[5];
  double rms[5];
  double mcmcnorm = 1./(double)NMCMCPts;
  for (int hbin=0;hbin<=hPassFail[0]->GetNbinsX();hbin++){
    mean[hbin]=0.;
    rms[hbin]=0.;
  }
  //calculate mean
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    for (int hbin=0;hbin<=hPassFail[0]->GetNbinsX();hbin++){
      mean[hbin] += mcmcnorm*hPassFail[ipt]->GetBinContent(hbin);
    }
  }
  mcmcnorm = 1./((double)NMCMCPts-1.);
  //calculate rms
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    for (int hbin=0;hbin<=hPassFail[0]->GetNbinsX();hbin++){
      double binc = hPassFail[ipt]->GetBinContent(hbin);
      double diff = binc-mean[hbin];
      rms[hbin] += TMath::Sqrt(diff*diff)*mcmcnorm;
    }
  }

  ////////////////////////////////////////////////////////
  //fill histogram with values
  for (int hbin=0;hbin<=hPassFail[0]->GetNbinsX();hbin++){
  //  cout<<"mean: "<<mean[hbin]<<endl;
  //  cout<<"RMS: "<<rms[hbin]<<endl;
  //  double binc = mean[hbin];
  //  double binerr = rms[hbin];
    hPassFailAvg->SetBinContent(hbin,mean[hbin]);
    hPassFailAvg->SetBinError(hbin,rms[hbin]);
  }


  //Draw everything on the same canvas
   hPassFailAvg->SetLineColor(kCyan+1);
//   hPassFailAvg->SetFillColor(kCyan+1);
   hPassFailAvg->SetLineWidth(3);
   hPassFailData->SetLineWidth(3);
   hPassFailData->SetMarkerStyle(8);
   hPassFail[0]->SetLineColor(kRed); //< this one is nominal MC
   hPassFail[0]->SetLineWidth(3);
   hPassFailData->Draw();
   hPassFail[0]->Draw("same");
   hPassFailAvg->Draw("same");
  
  return;

}

//void postfitCalculator::runPostFit(){

 // return;
//}


////////////////////////////////////////////////////////
//construct from parameters in parameter file
postfitCalculator::postfitCalculator(const char* parfile){
  
  //setup canvas
  TCanvas* cc = new TCanvas("cc","cc",800,700);

  //read in parameters
  runPars = new sharedPars(parfile);
  runPars->readParsFromFile();

  NMCMCPts = runPars->NMCMCPts;  
  MCMCBurnIn = runPars->MCMCBurnIn;
  NMCEvents = runPars->NMCEvents;
  NDataEvents = runPars->NDataEvents;

  //set parameter file name
  parFileName = parfile;

  //setup histogram manager
  hManager = new histoManager(runPars->hFactoryOutput.Data(),
                              runPars->nSamples,
                              runPars->nFVBins,
                              runPars->nComponents,
                              runPars->nAttributes);
  initHistos(); //< initilize attribute histograms using hmanager as a template 


  //get mcmc path tree
  TString mcmcfilename = runPars->MCMCFile;
  cout<<"postfitCalculator: getting mcmcpath from file: "<<mcmcfilename.Data()<<endl;
  TFile *mcmcfile = new TFile(mcmcfilename.Data());
  mcmcpath = (TTree*)mcmcfile->Get("MCMCpath");
  path = new mcmcReader(mcmcpath);

  //get mc and data trees
  cout<<"postfitCalculator: setting up data and MC trees from files:"<<endl;
  cout<<"    "<<runPars->hFactoryMCFiles.Data()<<endl;
  cout<<"    "<<runPars->hFactoryDataFiles.Data()<<endl;
  TChain* chmc = new  TChain("h1");
  chmc->Add(runPars->hFactoryMCFiles.Data());
  cout<<"postfitCalculator: # of MC events: "<<chmc->GetEntries()<<endl;
  TChain* chdat = new TChain("h1");
  chdat->Add(runPars->hFactoryDataFiles.Data());
  cout<<"postfitCalculator: # of Data events: "<<chdat->GetEntries()<<endl;;
  TTree* trmc = (TTree*)(chmc);
  TTree* trdata = (TTree*)(chdat);
  setMCTree(trmc);
  setDataTree(trdata);
  normFactor = (double)datatree->GetEntries()/(double)mctree->GetEntries(); 
  cout<<"postfitCalculator: normalization: "<<normFactor<<endl;
  //setup atmFitpars
  fitpars = new atmFitPars(parfile); 

  //get number of attributes in post fit
  nAttributes = runPars->nAttributes;  

  //initialize histograms
  init();    


}

/////////////////////////////////////////
//draws all modified histograms on same pad
void postfitCalculator::showAllHistos(int ihisto, int iclass){
  
  //loop over all histograms
  hEnuMu[iclass][0]->SetLineColor(kRed);
  if (ihisto==0)  hEnuMu[iclass][0]->Draw();
  if (ihisto==1)  hPIDemu[iclass][0]->Draw();
  for (int i=1;i<NMCMCPts;i++){
    if (ihisto==0){
      hEnuMu[iclass][i]->Draw("same");
    }
    if (ihisto==1){
      hPIDemu[iclass][i]->Draw("same");
    }
  }
  
  return;
}

/////////////////////////////////////////
//Fill table showing event numbers
void postfitCalculator::makeNSKTable(){

  //arrays to be filled
  double NSK[NCLASSMAX];
  double NSKRMS[NCLASSMAX];

  nsktree->SetBranchAddress("nevents",nevents);
  nsktree->SetBranchAddress("neventstot",neventstot);
  nsktree->SetBranchAddress("neventsdef",neventsdef);
  
  //set initial to zero
  for (int j=0;j<NCLASSMAX;j++){
    NSK[j]=0.;
    NSKRMS[j]=0.;
  }

  //loop over points to get mean
//  for (int i=0;i<nsktree->GetEntries();i++){
  for (int i=0;i<10;i++){
    nsktree->GetEntry(i);
    for (int j=0;j<NCLASSMAX;j++){
      NSK[j]+=nevents[j];
    }
  }
  for (int j=0;j<NCLASSMAX;j++){
    if (NSK[j]>=1.0){
      cout<<"NSK: "<<NSK[j]<<endl;
      NSK[j]/=10.; //< normalize for average # of events
    }
  }

  //loop over points to get RMS
//  for (int i=0;i<nsktree->GetEntries();i++){
  for (int i=0;i<10;i++){
    nsktree->GetEntry(i);
    for (int j=0;j<NCLASSMAX;j++){
      if (nevents[j]>0.1){
        NSKRMS[j] += ((nevents[j]-NSK[j])*(nevents[j]-NSK[j])); //< add square deviation
      }
    }
  }
  for (int j=0;j<NCLASSMAX;j++){
    if (NSKRMS[j]>=1.0){
      cout<<"NSKRMS: "<<NSKRMS[j]<<endl;      
      NSKRMS[j]/=10.; //< normalize 
      NSKRMS[j] = TMath::Sqrt( NSKRMS[j] ); //< get root mean square
    }
  }
  
  cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  for (int j=0;j<NCLASSMAX;j++){
    cout<<"Class: "<<j
    <<" Ndef: "<<neventsdef[j]
    <<" N: "<<NSK[j]
    <<" NRMS: "<<NSKRMS[j]
    <<" n: "<<(NSK[j]-neventsdef[j])/neventsdef[j]
    <<" nRMS: "<<NSKRMS[j]/neventsdef[j]
    <<endl;
  }
  cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;

  return;
}

///////////////////////////////////////////////////
//undo attribute modification
void postfitCalculator::unmodifyAttributes(){

  ////////////////////////////////////////////////////
  //un-modify attributes
  for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
    //mcreader->attribute[iatt]/=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0]; //divide by smear parameter
    mcreader->attribute[iatt]-=fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1]; //subtract bias parameter
  }

  //////////
  return;
}


///////////////////////////////////////////////////
//changes attributes of event to the modified quantities
void postfitCalculator::modifyAttributes(){
//  cout<<"N "<<fitpars->nAttributes<<endl;
  ////////////////////////////////////////////////////
  //modify attributes
  for (int iatt=0;iatt<fitpars->nAttributes;iatt++){
    int evtbin = mcreader->nbin;
    int evtcomp =  mcreader->ncomponent;
    int evtsamp = mcreader->nsample;
    double mean = hManager->hMCMean[evtsamp][evtbin][evtcomp][iatt];
    double scaling = fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][0];
    double bias    = fitpars->histoPar[mcreader->nbin][mcreader->ncomponent][iatt][1];
    mcreader->attribute[iatt] *= scaling; //multiply by smear parameter
    mcreader->attribute[iatt] += (mean - (scaling*mean)); //< correct mean shift from scaling   
    mcreader->attribute[iatt] += bias; //add bias parameter
  }

  //////////
  return;
}

///////////////////////////////////////////////
// modfies the values of the current MC event using parameters from current MCMC point
void postfitCalculator::modifyCurrentEvent(){
//  cout<<"test"<<endl;
  setParsFromMCMC(currentMCMCPoint); //< sets parameters from current step of MCMC path

  mctree->GetEntry(currentMCEvent); //< read in default MC event

  getEvtWeight(); //< gets the event weight form the current flux and xsec parameters

  modifyAttributes(); //< modifies each attribute according to the current MCMC point 
  
  return;
}

/////////////////////////////////////////////////////////////
//template for looping over MCMC points and modifying the MC
void postfitCalculator::MCMCLooperTemplate(){

  //get number of total steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-MCMCBurnIn;
  if (nsteps<(NMCMCPts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  //loop over points (uniform sampling)
  int dstep = nsteps/NMCMCPts;
  int istep = MCMCBurnIn;
  int samppts[10000];

  //get array of points to sample from
  for (int i=0;i<NMCMCPts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

  //find the number of MC events to use
  if ((NMCEvents<=0)||(NMCEvents>mctree->GetEntries())) NMCEvents = mctree->GetEntries(); 

  //loop over number of MCMC points
  for (int ipt=0;ipt<NMCMCPts;ipt++){
    currentMCMCPoint = samppts[ipt];
    //loop over MC events
    for (int ievt=0;ievt<NMCEvents;ievt++){
      currentMCEvent = ievt;
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      modifyCurrentEvent();//< all attributes are modified and eventWeight is calculated

      //feel free to do things here, like fill histograms or something  
 
    } 
  }
 
  return;

}



/////////////////////////////////////////////////
//Fill the NSK tree at various steps in MCMC path
void postfitCalculator::fillNskTree(int npts,int nburn){

  NMCMCPts = npts;
  
  //get number of steps minus burn-in
  int nsteps=mcmcpath->GetEntries()-nburn;
  if (nsteps<(npts)){
    cout<<"Not enough steps in MCMC path!"<<endl;
    return;
  }

  //loop over points (uniform sampling)
  int dstep = nsteps/npts;
  int istep = nburn;
  int samppts[10000];
  //get array of sampled points
  for (int i=0;i<npts;i++){
    samppts[i] = istep;
    istep+=dstep;
  } 

  //get total numbers of events
  findEvtClasses();


  //loop over number of points
  for (int ipt=0;ipt<npts;ipt++){
    setParsFromMCMC(samppts[ipt]);//< fill atmFitPars from parameters at this step
    for (int ievt=0;ievt<NMCEvents;ievt++){
      if ((ievt%1000)==0) cout<<"pt: "<<ipt<<" ev: "<<ievt<<endl;
      mctree->GetEntry(ievt);
      evtclass = getEvtClass(); //< get event class
      evtweight = getEvtWeight(); //< gets weight for event and modifies attributes   
      modifyAttributes(); //<  modify the attributes according the the mcmc parameters
      if (makeSelection(selectiontype)){
        fillHistos(evtclass,ipt); //< fills the histograms
        nevents[evtclass]++; //< increment the number of events of this class passing the cut
      }
    } 
 
    //fill nsk
    nsktree->Fill();

    //reset nevents
    for (int iclass=0;iclass<NCLASSMAX;iclass++){
      nevents[iclass]=0.;
    }
  }

  //write out output 
  nsktree->Write();
  //outFile->Close();

  ////////////////
  return; 
}

//////////////////////////////////////////////
//Find total numbers of events of each class
void postfitCalculator::findEvtClasses(){

  //first make sure each event class has zero events
  for (int iclass=0;iclass<NCLASSMAX;iclass++){
    neventstot[iclass]=0.;
    nevents[iclass]=0.;  
  }

  //now loop through all MC events to see hou many events of each
  //class we have
  cout<<"Finding total # of events of each class...."<<endl;
  for (int ievt=0;ievt<NMCEvents;ievt++){
    mctree->GetEntry(ievt);
    evtclass=getEvtClass(); //< fills "evtclass" variable
    neventstot[evtclass]++; //<increment the total number of events of this class
    if (makeSelection(selectiontype)){
      neventsdef[evtclass]++; //< increment the number of events of this class passing the cut
    }
  }
 
 
  /////////
  return;
}

//////////////////////////////////////////////
//Draw all event classes for a histogram
void postfitCalculator::drawBreakdown(int ihisto,int islot){

  ///////////////////
  //pid histograms
  if (ihisto==0){
    hPIDemu[0][islot]->SetFillColor(kCyan+1);
    hPIDemu[1][islot]->SetFillColor(kRed);
    hPIDemu[2][islot]->SetFillColor(kCyan);
    hPIDemu[3][islot]->SetFillColor(kOrange);
    hPIDemu[4][islot]->SetFillColor(kBlack); 
  }

  ///////////////////
  //Enu (mu) histograms
  if (ihisto==1){
    hEnuMu[0][islot]->SetFillColor(kCyan+1);
    hEnuMu[1][islot]->SetFillColor(kRed);
    hEnuMu[2][islot]->SetFillColor(kCyan);
    hEnuMu[3][islot]->SetFillColor(kOrange);
    hEnuMu[4][islot]->SetFillColor(kBlack);    
  }

  ////////////
  return;
}

/////////////////////////////////////////////////////////////
//returns an integer corresponding to the true event class
int postfitCalculator::getEvtClass(){
  if ((TMath::Abs(mcreader->mode)==1)&&(TMath::Abs(mcreader->ipnu[0])==12)) return 0; //< CCQE nu-e
  if ((TMath::Abs(mcreader->mode)==1)&&(TMath::Abs(mcreader->ipnu[0])==14)) return 1; //< CCQE nu-mu
  if ((TMath::Abs(mcreader->mode)>1)&&
      (TMath::Abs(mcreader->mode)<30)&&
      (TMath::Abs(mcreader->ipnu[0])==12)) return 2; //< CCnQE nu-e
  if ((TMath::Abs(mcreader->mode)>1)&&
      (TMath::Abs(mcreader->mode)<30)&&
      (TMath::Abs(mcreader->ipnu[0])==14)) return 3; //< CCnQE nu-mu
  if ((TMath::Abs(mcreader->mode)>30)) return 4; //< NC

  //////////////
  return -1;
}

void postfitCalculator::makeHistos(int iselect, int islot){
  
  //loop over mc events
  for (int ievt=0;ievt<mctree->GetEntries();ievt++){
    if ((ievt%1000)==0) cout<<ievt<<endl;
    mctree->GetEvent(ievt); //< load tree info
    evtweight = getEvtWeight(); //< gets weight for event and modifies attributes   
    modifyAttributes();
    evtclass = getEvtClass(); //< get event class
    if (makeSelection(iselect)){
      fillHistos(evtclass,islot); //< fills the histograms
      nevents[evtclass]++; //< increment the number of events of this class passing the cut
    }
  }

  //fill nsk tree for this point  
//  nsktree->Fill();

  ////////////
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//make various event selections, indexed by "iselect"
//returns 1 if all cuts are passed, 0 otherwise
int postfitCalculator::makeSelection(int iselect,int dataflg){

  fQreader* reader;
  if (dataflg) reader = datareader;
  else{
    reader = mcreader;   
  }

  ///////////////////////////
  //mock 1R numu calculation
  if (iselect==0){
    //Evis cut
    if (reader->attribute[2]<100.) return 0;
    //FC cut
    if (reader->nhitac>16) return 0;
    //PID cut
    if (reader->attribute[0]>0.) return 0;
    //nring cut
    if (reader->fqmrnring[0]!=1) return 0;
  }


  /////////////////////////////
  //cosmic mu mock selection
  if (iselect==1){
    //PID cut
    if ((-1.*reader->attribute[0])<0.2*reader->fq1rmom[0][1]) return 0;
    //nring cut
  //  if (mcreader->attribute[2]>0) return 0;
    //evis cut
    if (reader->fq1rmom[0][1]<30) return 0;
  }

  /////////////////////////////
  //cosmic decay e mock selection
  if (iselect==2){
    //PID cut
    if ((-1.*reader->attribute[1])>0.2*reader->fq1rmom[1][1]) return 0;
    //nring cut
  //  if (reader->attribute[2]>0) return 0;
    //evis cut
    if (reader->fq1rmom[0][1]<30) return 0;
  }

  //event has passed cuts
  return 1;
}


void postfitCalculator::fillHistos(int iclass,int islot){
  hEnuMu[iclass][islot]->Fill(mcreader->attribute[3],evtweight);
  hPIDemu[iclass][islot]->Fill(mcreader->attribute[0],evtweight);

}


///////////////////////////////////////////////////////////////////////////////////
//fills the atmFitPars object with the parameters from step "istep" of the MCMC
void  postfitCalculator::setParsFromMCMC(int istep){
  mcmcpath->GetEntry(istep); //< fills mcmcpath reader
  for (int ipar=0;ipar<path->npars;ipar++){
    fitpars->setParameter(ipar,path->par[ipar]);
  } 
  return;
}

void postfitCalculator::setDataTree(TTree* tr){
  datatree = tr;
  datareader = new fQreader(datatree);
  datatree->GetEntry(0);
  return;
}



void postfitCalculator::setMCTree(TTree* tr){
  mctree = tr;
  mcreader = new fQreader(mctree);
  mctree->GetEntry(0);
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//get event weight according to the current flux and xsec parameters
double  postfitCalculator::getEvtWeight(){

  /////////////////////////////////////////////
  //get event weights
  // double ww = mcreader->evtweight;
  double ww = 1.;
  //loop over all flux and xsec parameters
  for (int isyspar=0;isyspar<fitpars->nSysPars;isyspar++){
//  for (int isyspar=0;isyspar<1;isyspar++){
    double sysweight = sfact->getEvtWeight(mcreader,isyspar,fitpars->sysPar[isyspar]);
 //   cout<<"sysweight"<<sysweight<<endl; 
    ww*=sysweight;
  } 
  evtweight = ww;
  ///////////////
  return ww;
}


///////////////////////////////////////////////////////////////////////////////
//inialize some values
void postfitCalculator::init(){


  //////////////////////
  //make tree for numbers of events
  outFile = new TFile("nsk.root","recreate");
  nsktree = new TTree("nsk","nsk");
  nsktree->Branch("par",path->par,"par[100]/D");
  nsktree->Branch("nevents",nevents,"nevents[10]/D");
  nsktree->Branch("neventsdef",neventsdef,"neventsdef[10]/D");
  nsktree->Branch("neventstot",neventstot,"neventstot[10]/D");

  //////////////////////  
  //setup histograms
  int nbins=25;
  double xmin =  0.;
  double xmax = 5000.;
  double xminpid = -5000;
  double xmaxpid = 5000;
  for (int iclass=0;iclass<10;iclass++)
  {
    for (int i=0;i<200;i++){
      hEnuMu[iclass][i] = new TH1D(Form("enu_class%d_%d",iclass,i),Form("enu_class%d_%d",iclass,i) ,nbins,xmin,xmax);
      hPIDemu[iclass][i] = new TH1D(Form("PIDemu_class%d_%d",iclass,i),Form("PIDemu_class%d_%d",iclass,i) ,nbins,xminpid,xmaxpid);
    }
  }

  ////////////////////////////
  //setup atm fit pars
//  fitpars = new atmFitPars(parFileName.Data());
//  fitpars->initPars("tn186");

  ////////////////////////////
  //setup spline factory
  sfact = new splineFactory(parFileName.Data());  

  ////////////////////////////
  //other defaults
//  NMCEvents=0;
//  selectiontype=0;   
  return;
}


postfitCalculator::postfitCalculator(const char* mcmcfilename,const char* parfile){
  
  //set parameter file name
  parFileName = parfile;
  //get mcmc path
  TFile *mcmcfile = new TFile(mcmcfilename);
  mcmcpath = (TTree*)mcmcfile->Get("MCMCpath");
  path = new mcmcReader(mcmcpath);
  init();    
}

