#include "histoManager.h"

//TH1D*  histoManager::calcMCSum(int isample, int ibin, int iatt){
//  TH1D* hTot = (TH1D*)hMC[isample][ibin][0][iatt]->Clone("htot");
//  for (int i=1;i<nComponents;i++){
//    hTot->Add(hMC[isample][ibin][i][iatt]);
//  }
//  return hTot;
//}

///////////////////////////////////////////////////
//make some useful histograms for debugging
//------------do not use if separate neut mode
void histoManager::showSysParVariation(int isamp, int ibin, int icomp, int iatt, int ipar){

   //get array for parameter values to try
   const int npts = 15;
   double testval[npts];
   double minval =  fitPars->parUnc[ipar] - 3.*fitPars->parUnc[ipar];
   double maxval =  fitPars->parUnc[ipar] + 3.*fitPars->parUnc[ipar];
   double dval  = (maxval-minval)/((double)npts);
   double val = 0.;
   double currentval = fitPars->getParameter(ipar);
   for (int ipt=0;ipt<npts;ipt++){
     testval[ipt]= minval + val;
     val+=dval;
   } 
   cout<<"histoManager: Trying variations of param: "<<ipar<<" from "<<minval<<" to "<<maxval<<endl;

   //make 2D histogram
   int nbinsy = npts;
   TH1D* htemplate = getHistogram(isamp,ibin,icomp,iatt);
   int nbinsx = htemplate->GetNbinsX();
   double xbinmin = htemplate->GetBinLowEdge(1);
   double xbinmax = htemplate->GetBinLowEdge(nbinsx) + htemplate->GetBinWidth(nbinsx);
   h2d = new TH2D("h2d","h2d",nbinsx,xbinmin,xbinmax,nbinsy,minval,maxval);

   cout<<"filling histogram..."<<endl;
   //fill 2D histogram
   for (int ipt=0;ipt<npts;ipt++){
     fitPars->setParameter(ipar,testval[ipt]);
     cout<<"parval: "<<testval[ipt]<<endl;
     TH1D* htemp = getModHistogram(isamp,ibin,icomp,iatt);
     for (int b=1;b<=nbinsx;b++){
       h2d->SetBinContent(b,ipt,htemp->GetBinContent(b));
     } 
   }

   //draw 2D histogram
   h2d->GetZaxis()->SetRangeUser(0,h2d->GetMaximum());
   h2d->SetContour(50);
   h2d->Draw("lego2");
   fitPars->setParameter(ipar,currentval); 
   return;
}

//////////////////////////////////////////////////////////////////////////
//calculate sum of all modified components to compare to data
TH1D* histoManager::getSumHistogramMod(int isamp, int ibin, int iatt, int normFlg){
 
  /////////////////////////////////////////
  //erase previous calculation
  hSumHistoMod[isamp][ibin][iatt]->Reset();

  //////////////////////////////////////////
  //make sure this is set to properly account for errors
  hSumHistoMod[isamp][ibin][iatt]->SetDefaultSumw2(kTRUE);

  ///////////////////////////////////////////
  //add bin contents from all histograms
  if (!separateNeutMode) { // do not treat different interaction mode separately
    for (int icomp=0;icomp<nComponents;icomp++){
      TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,iatt);
      for (int jbin=1;jbin<=tmppointer->GetNbinsX();jbin++){
	double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
	//     content*=normFactor; //< scale the bin content
	content+=tmppointer->GetBinContent(jbin);
	double err1 = tmppointer->GetBinError(jbin);
	double err2 = hSumHistoMod[isamp][ibin][iatt]->GetBinError(jbin);
	hSumHistoMod[isamp][ibin][iatt]->SetBinContent(jbin,content);
	hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,TMath::Sqrt((err1*err1) + (err2*err2))); //<sum of squares weights
	// hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,100.); //<sum of squares weights
	//  hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,(err1+err2)); //<average scaling factors
      }        
    }
  } else { // treat different interaction mode separately
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int imode=0;imode<nModes;++imode) {
	TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,imode,iatt);
	for (int jbin=1;jbin<=tmppointer->GetNbinsX();jbin++){
	  double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
	  //     content*=normFactor; //< scale the bin content
	  content+=tmppointer->GetBinContent(jbin);
	  double err1 = tmppointer->GetBinError(jbin);
	  double err2 = hSumHistoMod[isamp][ibin][iatt]->GetBinError(jbin);
	  hSumHistoMod[isamp][ibin][iatt]->SetBinContent(jbin,content);
	  hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,TMath::Sqrt((err1*err1) + (err2*err2))); //<sum of squares weights
	} 
      }
    }
  }
  

  ////////////////////////////////////
  //return pointer to modified sum of componenets
  if (normFlg) hSumHistoMod[isamp][ibin][iatt]->Scale(normFactor); //< normalize the histo

  //calculate errors
//    for (int jbin=1;jbin<=hSumHistoMod[isamp][ibin][iatt]->GetNbinsX();jbin++){
//      double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
//      double err = hSumHistoMod[isamp][ibin][iatt]->GetBinError(jbin);
//      if (content>0.) hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,(err/content)); //<average scaling factors
//      else{
//        hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,0.);
//      }
//    }  

  return hSumHistoMod[isamp][ibin][iatt];

}

TH1D* histoManager::getSumHistogram(int isamp, int ibin, int iatt, int normFlg){
  hSumHisto[isamp][ibin][iatt]->Reset();
  hSumHisto[isamp][ibin][iatt]->SetDefaultSumw2(kTRUE);
  if (!separateNeutMode) {
    for (int icomp=0;icomp<nComponents;icomp++){
      hSumHisto[isamp][ibin][iatt]->Add(hMC[isamp][ibin][icomp][iatt]);
    }
  } else {
    for (int icomp=0;icomp<nComponents;icomp++) {
      for (int imode=0;imode<nModes; ++imode) {
	hSumHisto[isamp][ibin][iatt]->Add(hMCNeut[isamp][ibin][icomp][imode][iatt]);
      }
    }
  }

  if (normFlg) hSumHisto[isamp][ibin][iatt]->Scale(normFactor);
  return hSumHisto[isamp][ibin][iatt];
}

/////////////////////////////////////////////
//returns the given histogram after it has been
//modified by evaluating the spine weights
void histoManager::getSplineModifiedHisto(int isamp, int ibin, int icomp, int iatt){

  //number of bins in this thistogram
  int nhistobins = hMC[isamp][ibin][icomp][iatt]->GetNbinsX();

  //loop over bins and set new bin contents
  for (int i=1;i<=nhistobins;i++){
    double bincontent = hMC[isamp][ibin][icomp][iatt]->GetBinContent(i);
    double weightsum=0.;
    for (int isyspar=0;isyspar<fitPars->nSysPars;isyspar++){
      //get sum of spline weights from this bin
      weightsum *= getSplines(isamp,ibin,icomp,iatt)->evaluateSpline(i,isyspar,fitPars->sysPar[isyspar]);
    }
    //this formula gives the total weight to assign to this bin 
    //weightsum = weightsum -(double)fitPars->nSysPars + 1.; 
    bincontent*=weightsum;
    hMCModified[isamp][ibin][icomp][iatt]->SetBinContent(i,bincontent);
    binContents[i]=bincontent;
    hMCModified[isamp][ibin][icomp][iatt]->SetBinError(i,hMC[isamp][ibin][icomp][iatt]->GetBinError(i)*std::sqrt(weightsum));
  }
  
  //return  hMCModified[isamp][ibin][icomp][iatt];
}

void histoManager::getSplineModifiedHisto(int isamp, int ibin, int icomp, int imode, int iatt){
  //number of bins in this thistogram
  int nhistobins = hMCNeut[isamp][ibin][icomp][imode][iatt]->GetNbinsX();

  //loop over bins and set new bin contents
  for (int i=1;i<=nhistobins;i++){
    double bincontent = hMCNeut[isamp][ibin][icomp][imode][iatt]->GetBinContent(i);
    double weight = 1.0;
    for (int isyspar=0;isyspar<2;isyspar++){
      double wwgt = getSplines(isamp,ibin,icomp,imode,iatt)->evaluateSpline(i,isyspar,(fitPars->sysPar[isyspar]-fitPars->sysParNom[isyspar]));
      weight *= wwgt;//getSplines(isamp,ibin,icomp,imode,iatt)->evaluateSpline(i,isyspar,fitPars->sysPar[isyspar]);
      //std::cout<<fitPars->sysParName[isyspar]<<" "<<fitPars->sysPar[isyspar]<<" "<<(fitPars->sysPar[isyspar]-fitPars->sysParNom[isyspar])<<" "<<wwgt<<" | ";
    }
    for (int isysPar=3;isysPar<fitPars->nSysPars; ++isysPar) {
      if (fitPars->sysParName[isysPar].find("_C")!=std::string::npos) continue;
      if (fitPars->sysParName[isysPar].find("hc_")!=std::string::npos) continue;
      double wwgt;
      /*
      if (fitPars->sysParName[isysPar].find("RPA_O") != std::string::npos || 
	  fitPars->sysParName[isysPar].find("HAD") != std::string::npos ||
	  fitPars->sysParName[isysPar].find("MAQE") != std::string::npos ||
	  fitPars->sysParName[isysPar].find("pF_O") != std::string::npos ||
	  fitPars->sysParName[isysPar].find("EB_O") != std::string::npos || 
	  fitPars->sysParName[isysPar].find("CA5") != std::string::npos || 
	  fitPars->sysParName[isysPar].find("MANFFRES") != std::string::npos ||
	  fitPars->sysParName[isysPar].find("BgRES") != std::string::npos ||
	  fitPars->sysParName[isysPar].find("DISMPISHP") != std::string::npos ) {
	wwgt = getSplines(isamp,ibin,icomp,imode,iatt)->evaluateSpline(i,isysPar,(fitPars->sysPar[isysPar]-fitPars->sysParNom[isysPar]));
      } else {
	wwgt = getSplines(isamp,ibin,icomp,imode,iatt)->evaluateSpline(i,isysPar,fitPars->sysPar[isysPar]);
	}*/
      wwgt = getSplines(isamp,ibin,icomp,imode,iatt)->evaluateSpline(i,isysPar,(fitPars->sysPar[isysPar]-fitPars->sysParNom[isysPar])); 
      weight *= wwgt;//getSplines(isamp,ibin,icomp,imode,iatt)->evaluateSpline(i,isysPar,fitPars->sysPar[isysPar]);
      //std::cout<<fitPars->sysParName[isysPar]<<" "<<fitPars->sysPar[isysPar]<<" "<<(fitPars->sysPar[isysPar]-fitPars->sysParNom[isysPar])<<" "<<wwgt<<" | ";
    }
    bincontent*=weight;
    hMCNeutModified[isamp][ibin][icomp][imode][iatt]->SetBinContent(i,bincontent);
    binContents[i]=bincontent;
    hMCNeutModified[isamp][ibin][icomp][imode][iatt]->SetBinError(i,hMCNeut[isamp][ibin][icomp][imode][iatt]->GetBinError(i)*std::sqrt(weight));
  }
  
  //return  hMCNeutModified[isamp][ibin][icomp][imode][iatt];
}


TH1D* histoManager::getModHistogram(int isamp, int ibin, int icomp, int iatt){
  
  //number of bins in this histo
  int nhistobins = hMC[isamp][ibin][icomp][iatt]->GetNbinsX();
  double bincontent;

  //reset this histogram to erase any previous modifications
  hMCModified[isamp][ibin][icomp][iatt]->Reset();

  //apply splines if using them
  if (useSplineFlg){
    //hMCModified will now point to spline modified histogram
    getSplineModifiedHisto(isamp, ibin, icomp, iatt);
  }
  else{
    //just set bin contents equal to original histogram if not using splines
    for (int i=1;i<=nhistobins;i++){
      bincontent = hMC[isamp][ibin][icomp][iatt]->GetBinContent(i);
      binContents[i]=bincontent;
      hMCModified[isamp][ibin][icomp][iatt]->SetBinContent(i,bincontent);
      hMCModified[isamp][ibin][icomp][iatt]->SetBinError(i,hMC[isamp][ibin][icomp][iatt]->GetBinError(i));
    }
  }

  smearThisHistoFast( (*hMCModified[isamp][ibin][icomp][iatt]),binContents, fitPars->histoPar[ibin][icomp][iatt][0], fitPars->histoPar[ibin][icomp][iatt][1]);
  
  return hMCModified[isamp][ibin][icomp][iatt];
}

TH1D* histoManager::getModHistogram(int isamp, int ibin, int icomp, int imode, int iatt){ 
  //number of bins in this histo
  int nhistobins = hMCNeut[isamp][ibin][icomp][imode][iatt]->GetNbinsX();
  double bincontent;
  //reset this histogram to erase any previous modifications
  hMCNeutModified[isamp][ibin][icomp][imode][iatt]->Reset();
  //apply splines if using them
  if (useSplineFlg){
    //hMCModified will now point to spline modified histogram
    getSplineModifiedHisto(isamp, ibin, icomp, imode, iatt);
  }
  else{
    //just set bin contents equal to original histogram if not using splines
    for (int i=1;i<=nhistobins;i++){
      bincontent = hMCNeut[isamp][ibin][icomp][imode][iatt]->GetBinContent(i);
      binContents[i]=bincontent;
      hMCNeutModified[isamp][ibin][icomp][imode][iatt]->SetBinContent(i,bincontent);
      hMCNeutModified[isamp][ibin][icomp][imode][iatt]->SetBinError(i,hMCNeut[isamp][ibin][icomp][imode][iatt]->GetBinError(i));
    }
  }

  smearThisHistoFast( (*hMCNeutModified[isamp][ibin][icomp][imode][iatt]), binContents, fitPars->histoPar[ibin][icomp][iatt][0], fitPars->histoPar[ibin][icomp][iatt][1]);
  
  return hMCNeutModified[isamp][ibin][icomp][imode][iatt];
}


TH1D* histoManager::getHistogram(int isamp, int ibin, int icomp, int iatt){
  return hMC[isamp][ibin][icomp][iatt];
}

TH1D* histoManager::getHistogram(int isamp, int ibin, int icomp, int imode, int iatt){
  return hMCNeut[isamp][ibin][icomp][imode][iatt];
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
  if (separateNeutMode) {
    for (int i = 0; i < nComponents; ++i) {
      hMC[isample][ibin][i][iatt]->Reset();
      for (int j = 0; j < nModes; ++j) {
	for (int b = 1; b <= hMC[isample][ibin][i][iatt]->GetNbinsX(); ++b) 
	  hMC[isample][ibin][i][iatt]->SetBinContent(b, hMC[isample][ibin][i][iatt]->GetBinContent(b)+hMCNeut[isample][ibin][i][j][iatt]->GetBinContent(b));
      }
    }
  }
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
  if (!separateNeutMode) {
    for (int i = 0; i < nComponents; ++i) {
      hMC[isample][ibin][i][iatt]->Reset();
      for (int j = 0; j < nModes; ++j) {
	for (int b = 1; b <= hMC[isample][ibin][i][iatt]->GetNbinsX(); ++b) 
	  hMC[isample][ibin][i][iatt]->SetBinContent(b, hMC[isample][ibin][i][iatt]->GetBinContent(b)+hMCNeut[isample][ibin][i][j][iatt]->GetBinContent(b));
      }
    }
  }
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
histoManager::histoManager(const char* rootname,int nsamp,int nbin,int ncomp,int natt, int nmode, bool separateneutmode)
  : nModes(nmode)
  , separateNeutMode(separateneutmode)
{
  readFromFile(rootname,nsamp,nbin,ncomp,natt,nmode);
  nameTag = "histManager_For_";
  nameTag.Append(rootname);
  useSplineFlg=1;
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

void histoManager::setHistogram(int isamp, int ibin, int icomp, int imode, int iatt, int dataflg, TH1D *h)
{
  if (!dataflg){
    hMCNeut[isamp][ibin][icomp][imode][iatt] = h; 
  }
  else{
    hData[isamp][ibin][iatt] = h;
  }
}

void histoManager::fillHistogram(int isamp, int ibin, int icomp, int iatt, double value,double weight){
 // cout<<"sample: "<<isamp<<endl;
 // cout<<"bin"<<ibin<<endl;
 // cout<<"component"<<icomp<<endl;
 // cout<<"att"<<iatt<<endl;
 // cout<<"value"<<value<<endl;
 // cout<<"weight"<<weight<<endl;
  hMC[isamp][ibin][icomp][iatt]->Fill(value,weight); 
  return;
}

void histoManager::fillHistogram(int isamp, int ibin, int icomp, int imode, int iatt, double value, double weight)
{
  hMCNeut[isamp][ibin][icomp][imode][iatt]->Fill(value,weight);
}

void histoManager::fillHistogramData(int isamp, int ibin, int iatt, double value,double weight){
//  cout<<"sample: "<<isamp<<endl;
//  cout<<"bin"<<ibin<<endl;
//  cout<<"att"<<iatt<<endl;
//  cout<<"value"<<value<<endl;
//  cout<<"weight"<<weight<<endl;
  hData[isamp][ibin][iatt]->Fill(value,weight); 
  return;
}

//reads histograms from file
// ------- do not use if different neut modes are treated separately
void histoManager::readFromFile(const char* rootname,int nsamp,int nbin,int ncomp,int natt, int nmode)
{
  nSamples = nsamp;
  nBins    = nbin;
  nComponents = ncomp;
  nAttributes  = natt;
  //if (!nmode) nModes = nmode;
  nModes = NMODE;
  TString filename = rootname;

  //file containing histograms 
  fin = new TFile(filename.Data());
  ///////////////////////////////////////////////
  //get normalization factor between data and MC
  TH1D* htmp = (TH1D*)fin->Get("hnorm");

  ///////////////////////////////////////////////////
  //Set default sumw2 so errors are handled correctly
  //htmp->SetDefaultSumw2();
  TH1::SetDefaultSumw2(kTRUE);
  normFactor=htmp->GetBinContent(1);

  //////////////////////////////////////////////
  //read in histograms by name
  TString hname;
  if (!separateNeutMode) {
    //setup data histos
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
	for (int iatt=0;iatt<nAttributes;iatt++){
	  hname = "hdata_";
	  hname.Append(Form("samp%d_bin%d_att%d",isamp,ibin,iatt));
	  //cout<<"Getting histogram: "<<hname.Data()<<endl;
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
	    //cout<<"Getting histogram: "<<hname.Data()<<endl;
	    hMC[isamp][ibin][icomp][iatt] = (TH1D*)fin->Get(hname.Data());
	    //if (hMC[isamp][ibin][icomp][iatt]) cout<<"Got histogram: "<<hname.Data()<<endl;  
	    //        hmodname = hMC[isamp][ibin][icomp][iatt]->GetName();
	    //        hmodname.Append("_modified");
	    //        hMCModified[isamp][ibin][icomp][iatt] = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone(hmodname.Data());
	    //        hMCModified[isamp][ibin][icomp][iatt]->Reset();
	  }
	}
      }
    }
  } else { // treat different NEUT mode differently
    // data histos
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
	for (int iatt=0;iatt<nAttributes;iatt++){
	  hname = "hdata_";
	  hname.Append(Form("samp%d_bin%d_att%d",isamp,ibin,iatt));
	  //cout<<"Getting histogram: "<<hname.Data()<<endl;
	  hData[isamp][ibin][iatt] = (TH1D*)fin->Get(hname.Data());
	}
      }
    }
    //setup mc histos
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
	for (int icomp=0;icomp<nComponents;icomp++){
	  for (int imode = 0; imode < nModes; ++imode) {
	    for (int iatt=0;iatt<nAttributes;iatt++){
	      hname = "hmc_";
	      hname.Append(Form("samp%d_bin%d_comp%d_mode%d_att%d",isamp,ibin,icomp,imode,iatt));
	      //cout<<"Getting histogram: "<<hname.Data()<<endl;
	      hMCNeut[isamp][ibin][icomp][imode][iatt] = (TH1D*)fin->Get(hname.Data());
	    }
	  }
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
  TH1::SetDefaultSumw2(kTRUE);
  cout<<"histoManager: initializing modified histograms..."<<endl;
  TString hmodname;
  if (!separateNeutMode) {
    std::cout<<"do not separate neut mode"<<std::endl;
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
	for (int iatt=0;iatt<nAttributes;iatt++){
	  //if (hMC[isamp][ibin][0][iatt]==NULL) std::cout<<"what?"<<std::endl;
	  TString sumname = hMC[isamp][ibin][0][iatt]->GetName();
	  //std::cout<<sumname.Data()<<std::endl;
	  sumname.Append("_sum");
	  hSumHisto[isamp][ibin][iatt] = (TH1D*)hMC[isamp][ibin][0][iatt]->Clone(sumname.Data());
	  sumname.Append("_modified");
	  hSumHistoMod[isamp][ibin][iatt] = (TH1D*)hMC[isamp][ibin][0][iatt]->Clone(sumname.Data()); 
	  //std::cout<<sumname<<std::endl;
	  for (int icomp=0;icomp<nComponents;icomp++){
	    hmodname = hMC[isamp][ibin][icomp][iatt]->GetName();
	    hmodname.Append("_mod");
	    //std::cout<<hmodname<<std::endl;
	    hMCModified[isamp][ibin][icomp][iatt] = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone(hmodname.Data());
	    hMCModified[isamp][ibin][icomp][iatt]->Reset();
	  }
	}
      }
    }
  } else {
    std::cout<<"separate neut mode"<<std::endl;
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
	for (int iatt=0;iatt<nAttributes;iatt++){
	  TString sumname = hMCNeut[isamp][ibin][0][0][iatt]->GetName();
	  sumname.Append("_sum");
	  hSumHisto[isamp][ibin][iatt] = (TH1D*)hMCNeut[isamp][ibin][0][0][iatt]->Clone(sumname.Data());
	  sumname.Append("_modified");
	  hSumHistoMod[isamp][ibin][iatt] = (TH1D*)hMCNeut[isamp][ibin][0][0][iatt]->Clone(sumname.Data()); 
	  for (int icomp=0;icomp<nComponents;icomp++){
	    for (int imode = 0; imode < nModes; ++imode) {
	      hmodname = hMCNeut[isamp][ibin][icomp][imode][iatt]->GetName();
	      hmodname.Append("_mod");
	      //std::cout<<hmodname<<std::endl;
	      hMCNeutModified[isamp][ibin][icomp][imode][iatt] = (TH1D*)hMCNeut[isamp][ibin][icomp][imode][iatt]->Clone(hmodname.Data());
	      hMCNeutModified[isamp][ibin][icomp][imode][iatt]->Reset();
	    }
	  }
	}
      }
    }
  }
}

void histoManager::readSplinesFromFile(const char* fname,int nsyspartot){
  TString splinename;
  //make splines
  if (!separateNeutMode) {
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
  } else {
    for (int ibin=0;ibin<nBins;ibin++){
      for (int isamp=0;isamp<nSamples;isamp++){
	for (int icomp=0;icomp<nComponents;icomp++){
	  for (int imode = 0; imode < nModes; ++imode) {
	    for (int iatt=0;iatt<nAttributes;iatt++){
	      splinename="splinefor_";
	      splinename.Append(hMCNeut[isamp][ibin][icomp][imode][iatt]->GetName());
	      moreSplines[isamp][ibin][icomp][imode][iatt] = new hSplines(hMCNeut[isamp][ibin][icomp][imode][iatt],nsyspartot,splinename.Data());
	    }
	  }
	}
      }
    }
  }
  //open spline parameter file

  TFile splineFile(fname);
  TTree* splinePars = (TTree*)splineFile.Get("splinePars");
  splineParReader* parReader = new splineParReader(splinePars);

  //build the splines
  for (int ispline=0;ispline<splinePars->GetEntries();ispline++){
    splinePars->GetEntry(ispline);
    double Y[parReader->npoints];
    double X[parReader->npoints];
    for (int hbin=1;hbin<=parReader->nhistobins;hbin++){
      for (int kk=0;kk<parReader->npoints;kk++){
        Y[kk] = parReader->binWeight[kk][hbin];
        X[kk] = parReader->systParValues[kk];
	//std::cout<<X[kk]<<" "<<Y[kk]<<std::endl;
      }
      if (!separateNeutMode) {
	theSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nattribute]->buildSpline(hbin,parReader->nsystpar,X,Y,parReader->npoints);
      } else {
	//std::cout<<ispline<<" "<<hbin<<" "<<parReader->nsample<<" "<<parReader->nbin<<" "<<parReader->ncomponent<<" "<<parReader->nmode<<" "<<parReader->nattribute<<" "<<parReader->nsystpar<<" "<<parReader->npoints<<"-----------"<<std::endl;
        moreSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nmode][parReader->nattribute]->buildSpline(hbin, parReader->nsystpar, X, Y, parReader->npoints);
      }
    }
  }
  splineFile.Close();

  useSplineFlg=1;
  return;
}

histoManager::histoManager(int nsampl,int nbins,int ncomp,const char* name, int nmode, bool separateneutmode)
  : nModes(nmode)
  , separateNeutMode(separateneutmode)
{
  nameTag = "histos_";
  nameTag.Append(name);
  nSamples = nsampl;
  nComponents = ncomp;
  nAttributes = 0;
  nBins = nbins;
  useSplineFlg=1;
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
  nComponents = 2;
  nAttributes = 1;
  nBins = 1;
  useSplineFlg=1; 
  hMC[0][0][0][0]=testBumpD(nptsmc,3,0,"hmc_test");
  hMC[0][0][1][0]=testBumpD(nptsmc,3,0,"hmc_test_comp1");
  hMC[0][0][2][0]=testBumpD(nptsmc/3.,3,0,"hmc_test_comp2");
  hMC[0][0][3][0]=testBumpD(nptsmc,3,0,"hmc_test_comp3");
  hMC[0][0][4][0]=testBumpD(nptsmc,3,0,"hmc_test_comp4");
  hMC[0][0][5][0]=testBumpD(nptsmc,3,0,"hmc_test_comp5");

//  hMC[0][0][0][0]->SetName("mctest");
  hData[0][0][0]=testBumpD(nptsdata,3,0,"hdata_test");
  fitPars = new atmFitPars(1,1,nComponents,1,0);
  
  double nptstot = 0.;
  for (int i=0;i<nComponents;i++){
    nptstot += (double) hMC[0][0][i][0]->GetEntries();
  }
  fitPars->norm = (double)nptsdata/(nptstot);
  normFactor = fitPars->norm;
  hMC[0][0][0][0]->SetDefaultSumw2(kTRUE);
 // hMC[0][0][0][0]->Scale(fitPars->norm);
 // hMC[0][0][1][0]->Scale(fitPars->norm);
 // hMC[0][0][2][0]->Scale(fitPars->norm);
 // hMC[0][0][3][0]->Scale(fitPars->norm);
 // hMC[0][0][4][0]->Scale(fitPars->norm);
 // hMC[0][0][5][0]->Scale(fitPars->norm);

  initHistos();
  return;
}

void histoManager::showErrorComparison(int isamp, int ibin, int iatt){
   TH1D* hnom = getSumHistogram( isamp,  ibin,  iatt);
   TH1D* hmod = getSumHistogramMod( isamp, ibin, iatt);
   if (hTmp!=NULL) hTmp->Delete();
   hTmp = (TH1D*) hnom->Clone("herror");
   for (int jbin=1;jbin<=hnom->GetNbinsX();jbin++){
     double errdiff = TMath::Sqrt(hmod->GetBinContent(jbin)) - hmod->GetBinError(jbin);
     hTmp->SetBinContent(jbin,errdiff);
     hTmp->SetBinError(jbin,0.);
   }
  hTmp->Draw();
  return;
}

