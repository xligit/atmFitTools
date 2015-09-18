#ifndef HISTOMANAGER_C
#define HISTOMANAGER_C

#include "histoManager.h"


TH1F*  histoManager::calcMCSum(int isample, int ibin, int iatt){
  TH1F* hTot = (TH1F*)hMC[isample][ibin][0][iatt]->Clone("htot");
  for (int i=1;i<nComponents;i++){
    hTot->Add(hMC[isample][ibin][i][iatt]);
  }
  return hTot;
}



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
  float norm = 0.;
  if ((float)mcTree->GetEntries()>0) norm = (float)dataTree->GetEntries()/(float)mcTree->GetEntries();
  hMC[isample][ibin][hitolo[0]][iatt]->Draw();
  for (int j=0;j<nComponents;j++){
     hMC[isample][ibin][hitolo[j]][iatt]->Draw("same");
  }
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
  float norm = 0.;
  if ((float)mcTree->GetEntries()>0) norm = (float)dataTree->GetEntries()/(float)mcTree->GetEntries();
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

histoManager::histoManager(const char* rootname,int nsamp,int nbin,int ncomp,int natt){
  readFromFile(rootname,nsamp,nbin,ncomp,natt);
  nameTag = "histManager_For_";
  nameTag.Append(rootname);
  return;
}

void histoManager::readFromFile(const char* rootname,int nsamp,int nbin,int ncomp,int natt){
  nSamples = nsamp;
  nBins    = nbin;
  nComponents = ncomp;
  nAttributes  = natt;
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
//    if ((i%100)==0) cout<<"getting data event: "<<i<<endl;
    dataTree->GetEntry(i);
    fillAttributesData();
    for (int iatt=0;iatt<nAttributes;iatt++){
      hData[fqData->nsample][fqData->nbin][iatt]->Fill(att[iatt]);
   }
  }
  for (int j=0;j<nevmc;j++){
//    if ((j%100)==0) cout<<"getting mc event: "<<j<<endl;
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
  int nBinsNllEMu = 30;
  if (iatt==0){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-3000,6000);
  }
  if (iatt==1){
     hnew = new TH1F(thename,thename,nBinsNllEMu,-3000,6000);
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


#endif
