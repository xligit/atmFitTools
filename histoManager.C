#ifndef HISTOMANAGER_C
#define HISTOMANAGER_C

#include "histoManager.h"
#include "splineParReader.C"

//TH1D*  histoManager::calcMCSum(int isample, int ibin, int iatt){
//  TH1D* hTot = (TH1D*)hMC[isample][ibin][0][iatt]->Clone("htot");
//  for (int i=1;i<nComponents;i++){
//    hTot->Add(hMC[isample][ibin][i][iatt]);
//  }
//  return hTot;
//}



TH1D* histoManager::getSumHistogramMod(int isamp, int ibin, int iatt){
//  if (hSum){
  //  cout<<"delete previous sum histogram"<<endl;
 //   hSum->Delete();
 // }
//  hSum = getModHistogram(isamp,ibin,0,iatt);
  hSumHistoMod[isamp][ibin][iatt]->Reset();
  hSumHistoMod[isamp][ibin][iatt]->SetDefaultSumw2(kTRUE);
//  hSum->Sumw2(kTRUE);
//  hSum->Smooth(1);
 // cout<<"clone in new histogram"<<endl;
//  hSum = (TH1D*)getModHistogram(isamp,ibin,0,iatt)->Clone("hsum");
  for (int icomp=0;icomp<nComponents;icomp++){
    //cout<<"add histo component "<<icomp<<endl;
  //  tmppointer=getModHistogram(isamp,ibin,icomp,iatt);
  //  tmppointer->Smooth(3);
  //  hSum->Add(tmppointer);
//    tmppointer = getModHistogram(isamp,ibin,icomp,iatt);
//    for (int ibin=1;ibin<=hSum->GetNbinsX();ibin++){
//      double content = hSum->GetBinContent(ibin);
//      content+=tmppointer->GetBinContent(ibin);
//      double err1 = hSum->GetBinError(ibin);
//      double err2 = tmppointer->GetBinError(ibin);
//      hSum->SetBinContent(ibin,content);
 //     hSum->SetBinError(ibin,TMath::Sqrt((err1*err1) + (err2*err2)));
 //   }
 //   hSum->SetDefaultSumw2(kTRUE);
    hSumHistoMod[isamp][ibin][iatt]->Add(getModHistogram(isamp,ibin,icomp,iatt));
  }
  return hSumHistoMod[isamp][ibin][iatt];
}

TH1D* histoManager::getSumHistogram(int isamp, int ibin, int iatt){
//  if (hSumHisto[isamp][ibin][iatt]!=NULL) hSumHisto[isamp][ibin][iatt]->Delete();
  hSumHisto[isamp][ibin][iatt]->Reset();
  hSumHisto[isamp][ibin][iatt]->SetDefaultSumw2(kTRUE);
  for (int icomp=0;icomp<nComponents;icomp++){
    hSumHisto[isamp][ibin][iatt]->Add(hMC[isamp][ibin][icomp][iatt]);
  }
  return hSumHisto[isamp][ibin][iatt];
}

TH1D* histoManager::getModHistogram(int isamp, int ibin, int icomp, int iatt){
  int nhistobins = hMC[isamp][ibin][icomp][iatt]->GetNbinsX();
  double weightsum;
  double bincontent;
  hMCModified[isamp][ibin][icomp][iatt]->Reset();
  hMCModified[isamp][ibin][icomp][iatt]->SetDefaultSumw2(kTRUE);
  for (int i=1;i<=nhistobins;i++){
    bincontent = hMC[isamp][ibin][icomp][iatt]->GetBinContent(i);
    weightsum=0.;
    for (int isyspar=0;isyspar<fitPars->nSysPars;isyspar++){
      //debug
      //double daweight = getSplines(isamp,ibin,icomp,iatt)->evaluateSpline(i,isyspar,fitPars->sysPar[isyspar]);
      //cout<<daweight<<endl;
      //
      weightsum+=getSplines(isamp,ibin,icomp,iatt)->evaluateSpline(i,isyspar,fitPars->sysPar[isyspar]);
    } 
    weightsum = weightsum -(double)fitPars->nSysPars + 1.; 
    bincontent*=weightsum;
    hMCModified[isamp][ibin][icomp][iatt]->SetBinContent(i,bincontent);
    hMCModified[isamp][ibin][icomp][iatt]->SetBinError(i,hMC[isamp][ibin][icomp][iatt]->GetBinError(i));
  }
  smearThisHisto( (*hMCModified[isamp][ibin][icomp][iatt]), fitPars->histoPar[ibin][icomp][iatt][0], fitPars->histoPar[ibin][icomp][iatt][1]);
  return hMCModified[isamp][ibin][icomp][iatt];
}

/*
TH1D* histoManager::getModHistogram(int isamp, int ibin, int icomp, int iatt){
  if (hMod!=NULL){
 //   cout<<"deleting existing histogram"<<endl;
    hMod->Delete();
  }

  //if not using splines, start with base histogram
  if (!useSplineFlg){
  //  cout<<"getting base histogram"<<endl;
    hMod = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone("hmod");
  }

  //if splines are being used, start with histogram modified by splines
  else{
  //  cout<<"getting spline modified histogram"<<endl; 
    theSplines[isamp][ibin][icomp][iatt]->buildModHistoAllPar(fitPars->nSysPars,fitPars->sysPar);
//    theSplines[isamp][ibin][icomp][iatt]->buildModHisto(0,1.);

    hMod = (TH1D*)theSplines[isamp][ibin][icomp][iatt]->getModHisto()->Clone("hmod");
  }

  //apply histogram smearing
//  cout<<"smearing histogram with parameters: "<<fitPars->histoPar[ibin][0][iatt][0]<<" "<<fitPars->histoPar[ibin][0][iatt][1]<<endl;
//  TH1D* htmp = (TH1D*)hMod->CLone("htmo");
  smearThisHisto( (*hMod), fitPars->histoPar[ibin][icomp][iatt][0], fitPars->histoPar[ibin][icomp][iatt][1]);

  //return post-smearing histogram
  return hMod;
}
*/
TH1D* histoManager::getHistogram(int isamp, int ibin, int icomp, int iatt){
  return hMC[isamp][ibin][icomp][iatt];
}


//Use this to see how each MC component contributes to the overall histogram
void histoManager::showMCBreakdown(int isample,int ibin,int iatt){
  int color[NCOMPMAX];
  color[0] = 4;
  color[1] = 2;
  color[2] = 9;
  color[3] = 46;
  color[4] = 7;
  color[5] = 5;
  color[6] = 15;
  color[7] = 1;
  int style[NCOMPMAX];
  style[0] = 1001;
  style[1] = 1001;
  style[2] = 3001;
  style[3] = 3001;
  style[4] = 1001;
  style[5] = 1001;
  style[6] = 1001;
  style[7] = 1001;
  double size[NCOMPMAX];
  int hitolo[NCOMPMAX];
  for (int i=0;i<nComponents;i++){
//    hMC[isample][ibin][i][iatt]->SetLineColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillStyle(style[i]);
    size[i] = hMC[isample][ibin][i][iatt]->Integral();
    hitolo[i]=i;
  }
  int nswitch;
  //slow and easy 
  while (nswitch>0){
    nswitch=0;
    for (int ii=0;ii<(nComponents-1);ii++){
      if (size[hitolo[ii]]<size[hitolo[ii+1]]){
        nswitch = hitolo[ii];
        hitolo[ii] = hitolo[ii+1];
        hitolo[ii+1] = nswitch;
        nswitch=1;
      }
    }
  }
  double norm = 1.;
  //if ((double)mcTree->GetEntries()>0) norm = (double)dataTree->GetEntries()/(double)mcTree->GetEntries();
  hMC[isample][ibin][hitolo[0]][iatt]->Draw();
  for (int j=0;j<nComponents;j++){
     hMC[isample][ibin][hitolo[j]][iatt]->Draw("same");
  }
  if (Leg) Leg->Delete();
  Leg = new TLegend(0.7,0.6,0.9,0.9);
  Leg->AddEntry(hMC[isample][ibin][0][iatt],"CC1e","F");
  Leg->AddEntry(hMC[isample][ibin][1][iatt],"CC1#mu","F");
  Leg->AddEntry(hMC[isample][ibin][2][iatt],"CCeOth","F");
  Leg->AddEntry(hMC[isample][ibin][3][iatt],"CC#muOth","F");
  Leg->AddEntry(hMC[isample][ibin][4][iatt],"Single #pi^{0}","F");
  Leg->AddEntry(hMC[isample][ibin][5][iatt],"Single #pi^{+}","F");
  Leg->AddEntry(hMC[isample][ibin][6][iatt],"Other","F");
  Leg->Draw("same");
  return;
}

//stacks all MC histogram components to compare with data
THStack* histoManager::showMCBreakdownStack(int isample,int ibin,int iatt){
  int color[NCOMPMAX];
  color[0] = 4;
  color[1] = 2;
  color[2] = 9;
  color[3] = 46;
  color[4] = 7;
  color[5] = 5;
  color[6] = 15;
  color[7] = 1;
  int style[NCOMPMAX];
  style[0] = 1001;
  style[1] = 1001;
  style[2] = 3001;
  style[3] = 3001;
  style[4] = 1001;
  style[5] = 1001;
  style[6] = 1001;
  style[7] = 1001;
  double size[NCOMPMAX];
  int hitolo[NCOMPMAX];
  for (int i=0;i<nComponents;i++){
 //   hMC[isample][ibin][i][iatt]->SetLineColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillStyle(style[i]);
    size[i] = hMC[isample][ibin][i][iatt]->Integral();
    hitolo[i]=i;
  }
  int nswitch;
  //slow and easy 
  while (nswitch>0){
    nswitch=0;
    for (int ii=0;ii<(nComponents-1);ii++){
      if (size[hitolo[ii]]<size[hitolo[ii+1]]){
        nswitch = hitolo[ii];
        hitolo[ii] = hitolo[ii+1];
        hitolo[ii+1] = nswitch;
        nswitch=1;
      }
    }
  }
  double norm = 1.;
//  if ((double)mcTree->GetEntries()>0) norm = (double)dataTree->GetEntries()/(double)mcTree->GetEntries();
  THStack* hstack = new THStack("hstack","stack");
  hstack->Add(hMC[isample][ibin][hitolo[0]][iatt]);
  for (int j=1;j<nComponents;j++){
     hstack->Add(hMC[isample][ibin][hitolo[j]][iatt]);
  }
  hstack->Draw();
  hData[isample][ibin][iatt]->Scale(1./norm);
  hData[isample][ibin][iatt]->SetMarkerStyle(8);
  hData[isample][ibin][iatt]->Draw("samee");
  hData[isample][ibin][iatt]->Scale(norm);
  Leg = new TLegend(0.7,0.6,0.9,0.9);
  Leg->AddEntry(hMC[isample][ibin][0][iatt],"CC1e","F");
  Leg->AddEntry(hMC[isample][ibin][1][iatt],"CC1#mu","F");
  Leg->AddEntry(hMC[isample][ibin][2][iatt],"CCeOth","F");
  Leg->AddEntry(hMC[isample][ibin][3][iatt],"CC#muOth","F");
  Leg->AddEntry(hMC[isample][ibin][4][iatt],"Single #pi^{0}","F");
  Leg->AddEntry(hMC[isample][ibin][5][iatt],"Single #pi^{+}","F");
  Leg->AddEntry(hMC[isample][ibin][6][iatt],"Other","F");
  Leg->AddEntry(hData[isample][ibin][iatt],"Data","P");
  Leg->Draw("same");
  return hstack;
}

//constructor to re-created a histogram manager from a file
histoManager::histoManager(const char* rootname,int nsamp,int nbin,int ncomp,int natt){
  readFromFile(rootname,nsamp,nbin,ncomp,natt);
  nameTag = "histManager_For_";
  nameTag.Append(rootname);
  useSplineFlg=0;

  return;
}

void histoManager::setHistogram(int isamp, int ibin, int icomp, int iatt, int dataflg, TH1D* h){
  if (!dataflg){
    hMC[isamp][ibin][icomp][iatt] = h; 
  }
  else{
    hData[isamp][ibin][iatt] = h;
  }
  return;
}

void histoManager::fillHistogram(int isamp, int ibin, int icomp, int iatt, double value,double weight){
  hMC[isamp][ibin][icomp][iatt]->Fill(value,weight); 
  return;
}

void histoManager::fillHistogramData(int isamp, int ibin, int iatt, double value,double weight){
  hData[isamp][ibin][iatt]->Fill(value,weight); 
  return;
}

//reads histograms from file
void histoManager::readFromFile(const char* rootname,int nsamp,int nbin,int ncomp,int natt){
  nSamples = nsamp;
  nBins    = nbin;
  nComponents = ncomp;
  nAttributes  = natt;
  TString filename = rootname;
  //file containing histograms 
  fin = new TFile(filename.Data());
  ///////////////////////////////////////////////
  //get normalization factor between data and MC
  TH1D* htmp = (TH1D*)fin->Get("hnorm");

  ///////////////////////////////////////////////////
  //Set default sumw2 so errors are handled correctly
  htmp->SetDefaultSumw2();

  normFactor=htmp->GetBinContent(1);
  TString hname;
  //setup data histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int iatt=0;iatt<nAttributes;iatt++){
         hname = "hdata_";
         hname.Append(Form("samp%d_bin%d_att%d",isamp,ibin,iatt));
         cout<<"Getting histogram: "<<hname.Data()<<endl;
         hData[isamp][ibin][iatt] = (TH1D*)fin->Get(hname.Data());
      }
    }
  }
///  TString hmodname; //name for modified histogram
  //setup mc histos
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int icomp=0;icomp<nComponents;icomp++){
        for (int iatt=0;iatt<nAttributes;iatt++){
           hname = "hmc_";
           hname.Append(Form("samp%d_bin%d_comp%d_att%d",isamp,ibin,icomp,iatt));
           cout<<"Getting histogram: "<<hname.Data()<<endl;
           hMC[isamp][ibin][icomp][iatt] = (TH1D*)fin->Get(hname.Data());
   //        hmodname = hMC[isamp][ibin][icomp][iatt]->GetName();
   //        hmodname.Append("_modified");
   //        hMCModified[isamp][ibin][icomp][iatt] = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone(hmodname.Data());
   //        hMCModified[isamp][ibin][icomp][iatt]->Reset();
        }
      }
    }
  }
  initHistos(); 
  return;
}

///////////////////////////////////////////
//initializes arrays of modified histograms
void  histoManager::initHistos(){
  TH1D* hsumset;
  hsumset->SetDefaultSumw2(kTRUE);
  cout<<"histoManager: initializing modified histograms..."<<endl;
  TString hmodname;
  for (int isamp=0;isamp<nSamples;isamp++){
    for (int ibin=0;ibin<nBins;ibin++){
      for (int iatt=0;iatt<nAttributes;iatt++){
        TString sumname = hMC[isamp][ibin][0][iatt]->GetName();
        sumname.Append("_sum");
        hSumHisto[isamp][ibin][iatt] = (TH1D*)hMC[isamp][ibin][0][iatt]->Clone(sumname.Data());
        sumname.Append("_modified");
        hSumHistoMod[isamp][ibin][iatt] = (TH1D*)hMC[isamp][ibin][0][iatt]->Clone(sumname.Data()); 
        for (int icomp=0;icomp<nComponents;icomp++){
          hmodname = hMC[isamp][ibin][icomp][iatt]->GetName();
          hmodname.Append("_mod");
          hMCModified[isamp][ibin][icomp][iatt] = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone(hmodname.Data());
          hMCModified[isamp][ibin][icomp][iatt]->Reset();
        }
      }
    }
  }
  return; 
}

void histoManager::readSplinesFromFile(const char* fname,int nsyspartot){
  TString splinename;
  //make splines
  for (int ibin=0;ibin<nBins;ibin++){
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int icomp=0;icomp<nComponents;icomp++){
        for (int iatt=0;iatt<nAttributes;iatt++){
          splinename="splinefor_";
          splinename.Append(hMC[isamp][ibin][icomp][iatt]->GetName());
          theSplines[isamp][ibin][icomp][iatt]=new hSplines(hMC[isamp][ibin][icomp][iatt],nsyspartot,splinename.Data());
        }
      }
    }
  }
  //open spline parameter file
  TFile splineFile(fname);
  TTree* splinePars = (TTree*)splineFile.Get("splinePars");
  splineParReader* parReader = new splineParReader(splinePars);
  splinePars->GetEntry(0);
  //build the splines
  double Y[parReader->npoints];
  double X[parReader->npoints];
  for (int ispline=0;ispline<splinePars->GetEntries();ispline++){
    splinePars->GetEntry(ispline);
    for (int hbin=0;hbin<=parReader->nhistobins;hbin++){
      for (int ipt=0;ipt<parReader->npoints;ipt++){
        Y[ipt] = parReader->binWeight[ipt][hbin];
        X[ipt] = parReader->systParValues[ipt];
      }
//      cout<<parReader->nsample<<endl;
//      cout<<parReader->nbin<<endl;
//      cout<<parReader->ncomponent<<endl;
//      cout<<parReader->nattribute<<endl;

//      cout<< theSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nattribute]<<endl;
      theSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nattribute]
                  ->buildSpline(hbin,parReader->nsystpar,X,Y,parReader->npoints);

    }
  }
  splineFile.Close();
//  */
  useSplineFlg=1;
  return;
}

histoManager::histoManager(int nsampl,int nbins,int ncomp,const char* name){
  nameTag = "histos_";
  nameTag.Append(name);
  nSamples = nsampl;
  nComponents = ncomp;
  nAttributes = 0;
  nBins = nbins;
  useSplineFlg=0;
  return;
}

///////////////////////////////////////////////
//Unit testing
//Creates a histomanager object for gaussian histograms.
//'nptsmc' is the # of points in MC histos
//'nptsdata' is # of points in Data histos
//Useful for debugging likelihood evaluation
histoManager::histoManager(int nptsmc, int nptsdata){
  nameTag = "unittest";
  nSamples = 1;
  nComponents = 6;
  nAttributes = 1;
  nBins = 1;
  useSplineFlg=0; 
  hMC[0][0][0][0]=testBumpD(nptsmc,2,0,"hmc_test");
  hMC[0][0][1][0]=testBumpD(nptsmc,2,0,"hmc_test_comp1");
  hMC[0][0][2][0]=testBumpD(nptsmc,2,0,"hmc_test_comp2");
  hMC[0][0][3][0]=testBumpD(nptsmc,2,0,"hmc_test_comp3");
  hMC[0][0][4][0]=testBumpD(nptsmc,2,0,"hmc_test_comp4");
  hMC[0][0][5][0]=testBumpD(nptsmc,2,0,"hmc_test_comp5");

//  hMC[0][0][0][0]->SetName("mctest");
  hData[0][0][0]=testBumpD(nptsdata,2,0,"hdata_test");
  fitPars = new atmFitPars(1,1,nComponents,1,0);
  fitPars->norm = (double)nptsdata/((double)nptsmc*6.);
  hMC[0][0][0][0]->SetDefaultSumw2(kTRUE);
  hMC[0][0][0][0]->Scale(fitPars->norm);
  hMC[0][0][1][0]->Scale(fitPars->norm);
  hMC[0][0][2][0]->Scale(fitPars->norm);
  hMC[0][0][3][0]->Scale(fitPars->norm);
  hMC[0][0][4][0]->Scale(fitPars->norm);
  hMC[0][0][5][0]->Scale(fitPars->norm);

  initHistos();
  return;
}

void histoManager::showErrorComparison(int isamp, int ibin, int iatt){
   TH1D* hnom = getSumHistogram( isamp,  ibin,  iatt);
   TH1D* hmod = getSumHistogramMod( isamp, ibin, iatt);
   if (hTmp!=NULL) hTmp->Delete();
   hTmp = (TH1D*) hnom->Clone("herror");
   for (int ibin=1;ibin<=hnom->GetNbinsX();ibin++){
     double errdiff = TMath::Sqrt(hmod->GetBinContent(ibin)) - hmod->GetBinError(ibin);
     hTmp->SetBinContent(ibin,errdiff);
     hTmp->SetBinError(ibin,0.);
   }
  hTmp->Draw();
  return;
}

#endif
