#ifndef HISTOMANAGER_C
#define HISTOMANAGER_C

#include "histoManager.h"
#include "splineParReader.C"

//TH1F*  histoManager::calcMCSum(int isample, int ibin, int iatt){
//  TH1F* hTot = (TH1F*)hMC[isample][ibin][0][iatt]->Clone("htot");
//  for (int i=1;i<nComponents;i++){
//    hTot->Add(hMC[isample][ibin][i][iatt]);
//  }
//  return hTot;
//}



TH1F* histoManager::getSumHistogramMod(int isamp, int ibin, int iatt){
  if (hSum!=NULL) hSum->Delete();
  hSum = (TH1F*)getModHistogram(isamp,ibin,0,iatt)->Clone("hsum");
  for (int icomp=1;icomp<nComponents;icomp++){
    hSum->Add(getModHistogram(isamp,ibin,icomp,iatt));
  }
  return hSum;
}

TH1F* histoManager::getSumHistogram(int isamp, int ibin, int iatt){
  if (hSum!=NULL) hSum->Delete();
  hSum = (TH1F*)hMC[isamp][ibin][0][iatt]->Clone("hsum");
  for (int icomp=1;icomp<nComponents;icomp++){
    hSum->Add(hMC[isamp][ibin][icomp][iatt]);
  }
  return hSum;
}

TH1F* histoManager::getModHistogram(int isamp, int ibin, int icomp, int iatt){
  if (hMod!=NULL){
    cout<<"deleting existing histogram"<<endl;
    hMod->Delete();
  }

  //if not using splines, start with base histogram
  if (!useSplineFlg){
    cout<<"getting base histogram"<<endl;
    hMod = (TH1F*)hMC[isamp][ibin][icomp][iatt]->Clone("hmod");
  }

  //if splines are being used, start with histogram modified by splines
  else{
    cout<<"getting spline modified histogram"<<endl; 
    theSplines[isamp][ibin][icomp][iatt]->buildModHistoAllPar(fitPars->nSysPars,fitPars->sysPar);
//    theSplines[isamp][ibin][icomp][iatt]->buildModHisto(0,1.);

    hMod = (TH1F*)theSplines[isamp][ibin][icomp][iatt]->getModHisto()->Clone("hmod");
  }

  //apply histogram smearing
  cout<<"smearing histogram with parameters: "<<fitPars->histoPar[ibin][0][iatt][0]<<" "<<fitPars->histoPar[ibin][0][iatt][1]<<endl;
//  TH1F* htmp = (TH1F*)hMod->CLone("htmo");
//  smearThisHisto( (*hMod), fitPars->histoPar[ibin][0][iatt][0], fitPars->histoPar[ibin][0][iatt][1]);

  //return post-smearing histogram
  return hMod;
}

TH1F* histoManager::getHistogram(int isamp, int ibin, int icomp, int iatt){
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
  float size[NCOMPMAX];
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
  float norm = 1.;
  //if ((float)mcTree->GetEntries()>0) norm = (float)dataTree->GetEntries()/(float)mcTree->GetEntries();
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
  float size[NCOMPMAX];
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
  float norm = 1.;
//  if ((float)mcTree->GetEntries()>0) norm = (float)dataTree->GetEntries()/(float)mcTree->GetEntries();
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

void histoManager::setHistogram(int isamp, int ibin, int icomp, int iatt, int dataflg, TH1F* h){
  if (!dataflg){
    hMC[isamp][ibin][icomp][iatt] = h; 
  }
  else{
    hData[isamp][ibin][iatt] = h;
  }
  return;
}

void histoManager::fillHistogram(int isamp, int ibin, int icomp, int iatt, float value){
  hMC[isamp][ibin][icomp][iatt]->Fill(value); 
  return;
}

void histoManager::fillHistogramData(int isamp, int ibin, int iatt, float value){
  hData[isamp][ibin][iatt]->Fill(value); 
  return;
}

//reads histograms from file
void histoManager::readFromFile(const char* rootname,int nsamp,int nbin,int ncomp,int natt){
  nSamples = nsamp;
  nBins    = nbin;
  nComponents = ncomp;
  nAttributes  = natt;
  TString filename = rootname;
//  filename.Append(".root");
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

void histoManager::readSplinesFromFile(const char* fname){
  TFile splineFile(fname);
  TTree* splinePars = (TTree*)splineFile.Get("splinePars");
  splineParReader* parReader = new splineParReader(splinePars);
  splinePars->GetEntry(0);
  int nsyspartot = parReader->nsyspartot;
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
      theSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nattribute]
                  ->buildSpline(hbin,parReader->nsystpar,X,Y,parReader->npoints);
    }
  }
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


#endif
