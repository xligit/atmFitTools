#ifndef HISTOMANAGER_C
#define HISTOMANAGER_C

#include "histoManager.h"

//TH1D*  histoManager::calcMCSum(int isample, int ibin, int iatt){
//  TH1D* hTot = (TH1D*)hMC[isample][ibin][0][iatt]->Clone("htot");
//  for (int i=1;i<nComponents;i++){
//    hTot->Add(hMC[isample][ibin][i][iatt]);
//  }
//  return hTot;
//}

////////////////////////////////////////////////////
// Draw a particular spline
void histoManager::drawSpline(int isamp, int ibin, int icomp, int iatt, int hbin, int isyst){
  theSplines[isamp][ibin][icomp][iatt]->drawSpline(hbin,isyst);
  return;
}

////////////////////////////////////////////////////
// Draw 2D variation from splinesj
void histoManager::drawSpline2D(int isamp, int ibin, int icomp, int iatt, int isyst){
  int npts = 10;
  theSplines[isamp][ibin][icomp][iatt]->draw2D(npts,isyst);
  return;
}

///////////////////////////////////////////////////
//make some useful histograms for debugging
void histoManager::showSysParVariation(int isamp, int ibin, int icomp, int iatt, int ipar,double varscale){

   //get array for parameter values to try
   const int npts = 15;
   double testval[npts];
   double minval =  fitPars->pars[ipar] - 3.*varscale*fitPars->parUnc[ipar];
   double maxval =  fitPars->pars[ipar] + 3.*varscale*fitPars->parUnc[ipar];
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
     for (int jbin=1;jbin<=nbinsx;jbin++){
       h2d->SetBinContent(jbin,ipt,htemp->GetBinContent(jbin));
     } 
   }

   //draw 2D histogram
   if (h2d->GetMaximum()>0.) h2d->GetZaxis()->SetRangeUser(0,h2d->GetMaximum());
   h2d->SetContour(50);
   h2d->Draw("lego2");
   fitPars->setParameter(ipar,currentval); 
   return;
}

/*
//////////////////////////////////////////////////////////////////////////
//calculate sum of all modified components to compare to data
TH1D* histoManager::getSumHistogramNominal(int isamp, int ibin, int iatt, int normFlg){
 
  /////////////////////////////////////////
  //erase previous calculation
  hSumHistoMod[isamp][ibin][iatt]->Reset();

  //////////////////////////////////////////
  // set log likelihood to zero
  histoLogL=0.;

  //////////////////////////////////////////
  // number of edge bins to ignore
  int nbinbuffer = 3;
  nDOF = 0.;

  ///////////////////////////////////////////
  //add bin contents from all histograms
  if (!separateNeutMode) {
    for (int icomp=0;icomp<nComponents;icomp++){
      TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,iatt);
      for (int jbin=1;jbin<=tmppointer->GetNbinsX();jbin++){
      	double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
       	content+=tmppointer->GetBinContent(jbin);
        double normpar = fitPars->getNormParameter(isamp,ibin);
//        double normpar = 1.;
      	hSumHistoMod[isamp][ibin][iatt]->SetBinContent(jbin,normpar*content);
        if (icomp==(nComponents-1)){
        //  double normpar = fitPars->getNormParameter(isamp,ibin);
          if ((jbin>nbinbuffer) && jbin<(tmppointer->GetNbinsX()-nbinbuffer)){
            histoLogL+=evalLnL(hData[isamp][ibin][iatt]->GetBinContent(jbin),
                               normFactor*hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin));
            nDOF++;
          }
        }
      }        
    }

  }
  //histoLogL;

#ifdef T2K
  else {
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int imode=0;imode<nModes;++imode) {
        // get the modified histogram
      	TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,imode,iatt);
      	for (int jbin=1;jbin<=tmppointer->GetNbinsX();jbin++){
      	  double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
	        //     content*=normFactor; //< scale the bin content
       	  content+=tmppointer->GetBinContent(jbin);
      	  double err1 = tmppointer->GetBinError(jbin);
      	  double err2 = hSumHistoMod[isamp][ibin][iatt]->GetBinError(jbin);
      	  hSumHistoMod[isamp][ibin][iatt]->SetBinContent(jbin,content);
       	  hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,TMath::Sqrt((err1*err1) + (err2*err2))); //<sum of squares weights
      	  //std::cout<<"comp "<<icomp<<", mode "<<imode<<", jbin "<<jbin<<", "<<hSumHistoMod[isamp][ibin][iatt]->Integral()<<std::endl;
	      } 
      }
      // calculate likelihood after last component has been added
      if (icomp==(nComponents-1)){
        if ((jbin>nbinbuffer) && jbin<(tmppointer->GetNbinsX()-nbinbuffer)){
          histoLogL+=evalLnL(hData[isamp][ibin][iatt]->GetBinContent(jbin),
                               normFactor*hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin));
        }
      }
    }
  }
#endif

  ////////////////////////////////////
  //return pointer to modified sum of componenets
  if (normFlg) hSumHistoMod[isamp][ibin][iatt]->Scale(normFactor); //< normalize the histo


  return hSumHistoMod[isamp][ibin][iatt];

}

*/



//////////////////////////////////////////////////////////////////////////
//calculate sum of all modified components to compare to data
//NOTE: this also now calculates the log-likelihood
TH1D* histoManager::getSumHistogramMod(int isamp, int ibin, int iatt, int normFlg){
 
  /////////////////////////////////////////
  //erase previous calculation
  hSumHistoMod[isamp][ibin][iatt]->Reset();

  //////////////////////////////////////////
  // set log likelihood to zero
  histoLogL=0.;

  //////////////////////////////////////////
  // number of edge bins to ignore
  int nbinbuffer = 3;
  nDOF = 0.;

  ///////////////////////////////////////////
  //add bin contents from all histograms
  if (!separateNeutMode) {
    // loop over MC components
    for (int icomp=0;icomp<nComponents;icomp++){
      // pointer to temporary modified histogram
      TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,iatt);
      // loop over bins in this histogram
      for (int jbin=1;jbin<=tmppointer->GetNbinsX();jbin++){
        // get current summed content
      	double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
        // add temp histo content to sum
       	content+=tmppointer->GetBinContent(jbin);
        // properly normalize
        double normpar = 1.;
        if (useNormFlg) normpar = fitPars->getNormParameter(isamp,ibin);
        // set as new content
      	hSumHistoMod[isamp][ibin][iatt]->SetBinContent(jbin,normpar*content);
        // if on last component, compare to data and calc log L
        if (icomp==(nComponents-1)){
          if ((jbin>nbinbuffer) && jbin<(tmppointer->GetNbinsX()-nbinbuffer)){
            histoLogL+=evalLnL(hData[isamp][ibin][iatt]->GetBinContent(jbin),
                               normFactor*hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin));
//            cout<<"histoLogL: "<<histoLogL<<endl;
//            cout<<"samp: "<<isamp<<endl;
//            cout<<"bin: "<<ibin<<endl;
//            cout<<"att: "<<iatt<<endl;
            nDOF++;
          }
        }
      }        
    }

  }
  //histoLogL;

#ifdef T2K
  else {
    for (int icomp=0;icomp<nComponents;icomp++){
      for (int imode=0;imode<nModes;++imode) {
        // get the modified histogram
      	TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,imode,iatt);
      	for (int jbin=1;jbin<=tmppointer->GetNbinsX();jbin++){
      	  double content =  hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin);
	        //     content*=normFactor; //< scale the bin content
       	  content+=tmppointer->GetBinContent(jbin);
      	  double err1 = tmppointer->GetBinError(jbin);
      	  double err2 = hSumHistoMod[isamp][ibin][iatt]->GetBinError(jbin);
      	  hSumHistoMod[isamp][ibin][iatt]->SetBinContent(jbin,content);
       	  hSumHistoMod[isamp][ibin][iatt]->SetBinError(jbin,TMath::Sqrt((err1*err1) + (err2*err2))); //<sum of squares weights
      	  //std::cout<<"comp "<<icomp<<", mode "<<imode<<", jbin "<<jbin<<", "<<hSumHistoMod[isamp][ibin][iatt]->Integral()<<std::endl;
	      } 
      }
      // calculate likelihood after last component has been added
      if (icomp==(nComponents-1)){
        if ((jbin>nbinbuffer) && jbin<(tmppointer->GetNbinsX()-nbinbuffer)){
          histoLogL+=evalLnL(hData[isamp][ibin][iatt]->GetBinContent(jbin),
                               normFactor*hSumHistoMod[isamp][ibin][iatt]->GetBinContent(jbin));
        }
      }
    }
  }
#endif

  ////////////////////////////////////
  //return pointer to modified sum of componenets
  if (normFlg) hSumHistoMod[isamp][ibin][iatt]->Scale(normFactor); //< normalize the histo


  return hSumHistoMod[isamp][ibin][iatt];

}

////////////////////////////////////////////////////////////////////////////////////
// Returns NOMINAL summed histogram
TH1D* histoManager::getSumHistogram(int isamp, int ibin, int iatt, int normFlg){

    // reset sum histogram
    hSumHisto[isamp][ibin][iatt]->Reset();
    hSumHisto[isamp][ibin][iatt]->SetDefaultSumw2(kTRUE);
  
    // loop over components and add them
    for (int icomp=0;icomp<nComponents;icomp++){

      // get current parameters
      double tmpsmear = fitPars->getHistoParameter(ibin,icomp,iatt,0);
      double tmpbias = fitPars->getHistoParameter(ibin,icomp,iatt,1);
      int tmpsplineflg = useSplineFlg;
      int tmpnormflg   = useNormFlg; 

      // find which parameter effect this histo
      int smearindex = fitPars->getParIndex(ibin,icomp,iatt,0);
      int biasindex = fitPars->getParIndex(ibin,icomp,iatt,1);

      // reset these pars to nominal values
      fitPars->setParameter(smearindex,1.0);
      fitPars->setParameter(biasindex,0.0);
      useSplineFlg = 0;
      useNormFlg = 0;

      // get temporary pointer to histogra
      // (this does not create a new histogram)
      TH1D* tmppointer=getModHistogram(isamp,ibin,icomp,iatt);

      // add histogram to total
      hSumHisto[isamp][ibin][iatt]->Add(tmppointer);
      
      // restore current parameters
      fitPars->setParameter(smearindex,tmpsmear);
      fitPars->setParameter(biasindex,tmpbias);
      useSplineFlg = tmpsplineflg;
      useNormFlg = tmpnormflg;
  }

  if (normFlg) hSumHisto[isamp][ibin][iatt]->Scale(normFactor);
  return hSumHisto[isamp][ibin][iatt];

  /*
  hSumHisto[isamp][ibin][iatt]->Reset();
  hSumHisto[isamp][ibin][iatt]->SetDefaultSumw2(kTRUE);
  TH1D* htmp = (TH1D*)hSumHisto[isamp][ibin][iatt]->Clone("htmp_sumhisto");
  if (!separateNeutMode) {
    for (int icomp=0;icomp<nComponents;icomp++){
      rebinHisto( hMC[isamp][ibin][icomp][iatt],htmp);
      hSumHisto[isamp][ibin][iatt]->Add(htmp);

    }
  }
  htmp->Delete();
#ifdef T2K
  else {
    for (int icomp=0;icomp<nComponents;icomp++) {
      for (int imode=0;imode<nModes; ++imode) {
	hSumHisto[isamp][ibin][iatt]->Add(hMCNeut[isamp][ibin][icomp][imode][iatt]);
      }
    }
  }
#endif

  if (normFlg) hSumHisto[isamp][ibin][iatt]->Scale(normFactor);
  return hSumHisto[isamp][ibin][iatt];
  */
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the spline-modified contents of the specified histogram bin
double histoManager::getSplineModifiedBin(int isamp, int ibin, int icomp, int iatt, int ihistobin){
 
 if (ihistobin<1) return 0.;
 if (ihistobin>hMC[isamp][ibin][icomp][iatt]->GetNbinsX()) return 0.;

 // originial bin contents
 double bincontent = hMC[isamp][ibin][icomp][iatt]->GetBinContent(ihistobin);
// cout<<"content: "<<bincontent<<endl;
 // sum of weights from splines
 double weightsum=0.;
 for (int isyspar=0;isyspar<(fitPars->nSysPars-fitPars->nNormPars);isyspar++){
//    cout<<"weight: "<<weightsum<<endl;
    weightsum+= getSplines(isamp,ibin,icomp,iatt)->evaluateSpline(ihistobin,isyspar,fitPars->sysPar[isyspar]);
//    weightsum+= theSplines[isamp][ibin][icomp][iatt]->evaluateSpline(ihistobin,isyspar,fitPars->sysPar[isyspar]);
  }
  //this formula gives the total weight to assign to this bin 
  weightsum = weightsum -(double)fitPars->nSysPars + (double)fitPars->nNormPars + 1.; 
//  cout<<"weight final: "<<weightsum<<endl;
  if (weightsum>0.) bincontent*=weightsum;
 
  //
  return  bincontent;

}

/////////////////////////////////////////////
//returns the given histogram after it has been
//modified by evaluating the spine weights
TH1D* histoManager::getSplineModifiedHisto(int isamp, int ibin, int icomp, int iatt){


  //number of bins in this thistogram
  int nhistobins = hMC[isamp][ibin][icomp][iatt]->GetNbinsX();

  //loop over bins and set new bin contents
  for (int i=1;i<=nhistobins;i++){
    double bincontent = hMC[isamp][ibin][icomp][iatt]->GetBinContent(i);
    double weightsum=0.;
    for (int isyspar=0;isyspar<fitPars->nSysPars;isyspar++){
      //get sum of spline weights from this bin
//      cout<<"evaluating: spline for bin: "<<i<<endl;
 //     cout<<"for histo : "<< hMC[isamp][ibin][icomp][iatt]->GetName()<<endl;;
      weightsum+=getSplines(isamp,ibin,icomp,iatt)->evaluateSpline(i,isyspar,fitPars->sysPar[isyspar]);
    }
    //this formula gives the total weight to assign to this bin 
    weightsum = weightsum -(double)fitPars->nSysPars + (double)fitPars->nNormPars + 1.; 
    bincontent*=weightsum;
    hMCModified[isamp][ibin][icomp][iatt]->SetBinContent(i,bincontent);
    binContents[i]=bincontent;
    hMCModified[isamp][ibin][icomp][iatt]->SetBinError(i,hMC[isamp][ibin][icomp][iatt]->GetBinError(i));
  }
  
  return  hMCModified[isamp][ibin][icomp][iatt];
}

#ifdef T2K
TH1D* histoManager::getSplineModifiedHisto(int isamp, int ibin, int icomp, int imode, int iatt){
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
      //std::cout<<fitPars->sysParName[isysPar]<<" "<<fitPars->sysPar[isysPar]<<" "<<(fitPars->sysPar[isysPar]-fitPars->sysParNom[isysPar])<<" "<<wwgt<<"\n";
    }
    bincontent*=weight;
    hMCNeutModified[isamp][ibin][icomp][imode][iatt]->SetBinContent(i,bincontent);
    binContents[i]=bincontent;
    hMCNeutModified[isamp][ibin][icomp][imode][iatt]->SetBinError(i,hMCNeut[isamp][ibin][icomp][imode][iatt]->GetBinError(i)*TMath::Sqrt(weight));
  }
  
  return  hMCNeutModified[isamp][ibin][icomp][imode][iatt];
}

#endif

//////////////////////////////////////////////////////////////////////////////
// Fills gTmp a modifed GRAPH of the specified histogram
TGraph* histoManager::getModGraph(int isamp, int ibin, int icomp, int iatt){

  // goal integral of graph (sum of bins)
  gTotIntegral = 0.;

  // get modification parameters
  double bias = fitPars->getHistoParameter(ibin,icomp,iatt,1);
  double smear = fitPars->getHistoParameter(ibin,icomp,iatt,0);
  if (smear<0.) smear*=-1.0;

  // number of bins in this graph
  int nbins = hMC[isamp][ibin][icomp][iatt]->GetNbinsX();

  // make graph 
  TGraph* gr = new TGraph(nbins+2); 
 
  // fill graph
  double binc = 0; //< bin contens
  double binw = hMC[isamp][ibin][icomp][iatt]->GetBinWidth(1); //< bin width

  // shift from smear par
  //double maxx = hMCMean[isamp][ibin][icomp][iatt];
 // double smearshift = (maxx*smear) - maxx;
 
  gr->SetPoint(0,hMC[isamp][ibin][icomp][iatt]->GetBinLowEdge(1),0);
  for (int ipt=1; ipt<=nbins; ipt++){
    
    // get spline modified contents if necessary
    if (useSplineFlg){
      binc = getSplineModifiedBin(isamp,ibin,icomp,iatt,ipt);
    }
    else{
      binc = hMC[isamp][ibin][icomp][iatt]->GetBinContent(ipt);
    }

    // add to total integral
    gTotIntegral+= binc;

    // scale bin contents by scaling parameter to keep graph integral the same
    binc/=smear;

    // scale by bin width for normalization
    binc/=binw;

    // fill modified graph
    double xpt =  ((smear*hMC[isamp][ibin][icomp][iatt]->GetBinCenter(ipt)) + bias);
//    double xpt =  ((smear*hMC[isamp][ibin][icomp][iatt]->GetBinLowEdge(ipt)+binw) + bias);
    double ypt = binc;
    gr->SetPoint(ipt,xpt,ypt);

//    gr->SetPoint(ipt,
//                 ((smear*hMC[isamp][ibin][icomp][iatt]->GetBinCenter(ipt)) + bias),
//                 binc);
  }

  gr->SetPoint(nbins+1,
              (hMC[isamp][ibin][icomp][iatt]->GetBinLowEdge(nbins)+hMC[isamp][ibin][icomp][iatt]->GetBinWidth(nbins))*smear + bias,0);
 
  // smooth it?
//  smoothGraph(gr);

  return gr;
}


/////////////////////////////////////////////////////////////////////////////
// Apply a lower bound to some attribute
void histoManager::setLoBound(int iatt, double bound){

  // make changes
  applyPhysicalBound[iatt] = 1;
  physicalLoBound[iatt] = bound;
 
  // talk about it
  cout<<"histoManager: Adding lower bound of "<<bound<<"to parameter "<<iatt<<endl;
  return;
}

//////////////////////////////////////////////////////////////////////////////
// NEW method using graphs
// Returns the modified histogram based on the parameters in atmfitpars
TH1D* histoManager::getModHistogram(int isamp, int ibin, int icomp, int iatt){

  // get the modified graph and calculate gTotIntegral
  TGraph* gr = getModGraph(isamp, ibin, icomp, iatt);
 
  // get norm factor
  double histonorm = gTotIntegral;

  // convert to to histogram 
  graph2histo(gr,hMCModified[isamp][ibin][icomp][iatt]);

  // apply any physical bounds
  if (applyPhysicalBound[iatt]){
    applyLoBound(gr,hMCModified[isamp][ibin][icomp][iatt]);
  }

  
  gr->Delete();
  return hMCModified[isamp][ibin][icomp][iatt];
}


//////////////////////////////////////////////////////////////////////////////
// NEW method using graphs 
// Returns the modified histogram based on the parameters in atmfitpars
// This function return the mod histogram in original binning
TH1D* histoManager::getModHistogramMC(int isamp, int ibin, int icomp, int iatt){

  // get the modified graph and calculate gTotIntegral
  TGraph* gr = getModGraph(isamp, ibin, icomp, iatt);
 
  // convert to to histogram 
  TH1D* htmp = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone("htmp");
  graph2histo(gr,htmp);

  // apply physical bounds if any
  if (applyPhysicalBound[iatt]){
    applyLoBound(gr,htmp);
  }

  gr->Delete();
  return htmp;
}



/*
//////////////////////////////////////////////////////////////////////////////
// Returns the modified histogram based on the parameters in atmfitpars
TH1D* histoManager::getModHistogram(int isamp, int ibin, int icomp, int iatt){

  //number of bins in this histo
  int nhistobins = hMC[isamp][ibin][icomp][iatt]->GetNbinsX();

  //reset this histogram to erase any previous modifications
  hMCModified[isamp][ibin][icomp][iatt]->Reset();

  //bin width
  double binw = hMC[isamp][ibin][icomp][iatt]->GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double sumw;
  double binc;
  double bias = fitPars->getHistoParameter(ibin,icomp,iatt,1);
  double smear = 1./fitPars->getHistoParameter(ibin,icomp,iatt,0);
//  double norm = fitPars->getNormParameter(isamp,ibin)


  // smooth weights
  double ww0 = 3.98942e-01;
  double ww1 = 5.39909e-02;
  double ww2 = 1.33830e-04;
  double wwsum = ww0 + ww1 + ww2;
  double contsum=0.;
  double wcontsum=0.;

  // keep a list of bin contents
  deque<double> content_list;

  // initialize list
  if (!useSplineFlg){
    content_list.push_back(0.);
    content_list.push_back(0.);
    content_list.push_back(hMC[isamp][ibin][icomp][iatt]->GetBinContent(1));
    content_list.push_back(hMC[isamp][ibin][icomp][iatt]->GetBinContent(2));
    content_list.push_back(hMC[isamp][ibin][icomp][iatt]->GetBinContent(3));
  }
  else{
    content_list.push_back(0.);
    content_list.push_back(0.);
    content_list.push_back(getSplineModifiedBin(isamp,ibin,icomp,iatt,1));
    content_list.push_back(getSplineModifiedBin(isamp,ibin,icomp,iatt,2));
    content_list.push_back(getSplineModifiedBin(isamp,ibin,icomp,iatt,3));
  }
  //loop over bins and modify contents
  for (int newbin=1;newbin<=nhistobins;newbin++){
    sum = 0.;
    sumw = 0;
    binerr=0.;
    binedge = hMCModified[isamp][ibin][icomp][iatt]->GetBinLowEdge(newbin);
    ymin = (binedge - bias)*smear;
    ymax = ymin + binw;
    for (int oldbin=1;oldbin<=nhistobins;oldbin++){
      xmin = hMCModified[isamp][ibin][icomp][iatt]->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      // calculate smoothed bin contents
      binc =  ww0*content_list.at(2);
      contsum+=(binc/ww0);
      binc += ww1*content_list.at(3);
      binc += ww1*content_list.at(1);
      binc += ww2*content_list.at(4);
      binc += ww2*content_list.at(0);
      binc/=(wwsum);
      wcontsum+=binc;
      //update list
      content_list.pop_front();
      if ((oldbin+3)>nhistobins) content_list.push_back(0.);
      else if (useSplineFlg){
        content_list.push_back(getSplineModifiedBin(isamp,ibin,icomp,iatt,oldbin+3));
      }
      else{
        content_list.push_back(hMC[isamp][ibin][icomp][iatt]->GetBinContent(oldbin+3));
      }      
      sum+=weight*binc;
      binerr += weight*weight*binc;
    }
    hMCModified[isamp][ibin][icomp][iatt]->SetBinContent(newbin,sum);
    hMCModified[isamp][ibin][icomp][iatt]->SetBinError(newbin,TMath::Sqrt(binerr));
  }

  hMCModified[isamp][ibin][icomp][iatt]->Scale(fitPars->getNormParameter(isamp,ibin)*contsum/wcontsum);

  return hMCModified[isamp][ibin][icomp][iatt];
  
}
*/


//////////////////////////////////////////////////////////////////////////////
// Returns the modified histogram based on the parameters in atmfitpars
TH1D* histoManager::getModHistogramSlow(int isamp, int ibin, int icomp, int iatt){
  
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
  // do nothing if histogram has hardly any events
  if (hMC[isamp][ibin][icomp][iatt]->GetEntries()<10) return hMC[isamp][ibin][icomp][iatt];
  // otherwise, modifiy this histogram
//  smearThisHistoFast( (*hMCModified[isamp][ibin][icomp][iatt]),
//                        binContents, 
//                        fitPars->getHistoParameter(ibin,icomp,iatt,0),
//                        fitPars->getHistoParameter(ibin,icomp,iatt,1),
//			                  fitPars->getNormParameter(isamp,ibin) );
  // otherwise, modifiy this histogram with mean
  smearThisHistoFastMean( (*hMCModified[isamp][ibin][icomp][iatt]),
                        binContents, 
                        fitPars->getHistoParameter(ibin,icomp,iatt,0),
                        hMCMean[isamp][ibin][icomp][iatt],
                        fitPars->getHistoParameter(ibin,icomp,iatt,1),
			                  fitPars->getNormParameter(isamp,ibin) );
  
  return hMCModified[isamp][ibin][icomp][iatt];
}



#ifdef T2K
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
  //smearThisHistoFast( (*hMCNeutModified[isamp][ibin][icomp][imode][iatt]), binContents, fitPars->histoPar[ibin][icomp][iatt][1]);

  smearThisHistoFast( (*hMCNeutModified[isamp][ibin][icomp][imode][iatt]), binContents, fitPars->histoPar[ibin][icomp][iatt][0], fitPars->histoPar[ibin][icomp][iatt][1]);

  return hMCNeutModified[isamp][ibin][icomp][imode][iatt];
}

TH1D* histoManager::getHistogram(int isamp, int ibin, int icomp, int imode, int iatt){
  return hMCNeut[isamp][ibin][icomp][imode][iatt];
}

TH1D *histoManager::getNominalHistogram(int isamp, int ibin, int icomp, int imode, int iatt)
{
  return hMCNeutNom[isamp][ibin][icomp][imode][iatt];
}

#endif

TH1D* histoManager::getHistogram(int isamp, int ibin, int icomp, int iatt){
  return hMC[isamp][ibin][icomp][iatt];
}

/////////////////////////////////////////////////////////////////////////////
//Use this to see how each MC component contributes to the overall histogram
void histoManager::showMCBreakdown(int isample,int ibin,int iatt){

  // setup plots, colors, etc.
  int color[NCOMPMAX];
  color[0] = 4;
  color[1] = 2;
  color[2] = 7;
  color[3] = 6;
  color[4] = 5;
  color[5] = 8;
  //color[6] = 15;
//  color[7] = 1;
  int style[NCOMPMAX];
  style[0] = 1001;
  style[1] = 1001;
  style[2] = 1001;
  style[3] = 1001;
  style[4] = 1001;
  style[5] = 1001;
//  style[6] = 1001;
//  style[7] = 1001;
  double size[NCOMPMAX];
  int hitolo[NCOMPMAX];

  for (int i=0;i<nComponents;i++){
    hMC[isample][ibin][i][iatt]->SetFillColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillStyle(style[i]);
    size[i] = hMC[isample][ibin][i][iatt]->Integral();
    hitolo[i]=i;
  }

  int nswitch=1;
  //slow but easy 
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

  // now make plot
  hMC[isample][ibin][hitolo[0]][iatt]->Draw("h");
  for (int j=1;j<nComponents;j++){
     hMC[isample][ibin][hitolo[j]][iatt]->Draw("sameh");
  }
  if (Leg) Leg->Delete();
  Leg = new TLegend(0.7,0.6,0.9,0.9);
  Leg->AddEntry(hMC[isample][ibin][0][iatt],"Single e","F");
  Leg->AddEntry(hMC[isample][ibin][1][iatt],"Single #mu","F");
  Leg->AddEntry(hMC[isample][ibin][2][iatt],"e + Other","F");
  Leg->AddEntry(hMC[isample][ibin][3][iatt],"#mu + Other","F");
  Leg->AddEntry(hMC[isample][ibin][4][iatt],"Single #pi^{0}","F");
//  Leg->AddEntry(hMC[isample][ibin][5][iatt],"Other","F");
 // Leg->AddEntry(hMC[isample][ibin][6][iatt],"Other","F");
  Leg->Draw("same");
  return;
}

//stacks all MC histogram components to compare with data
THStack* histoManager::showMCBreakdownStack(int isample,int ibin,int iatt){
  int color[NCOMPMAX];
  color[0] = 4;
  color[1] = 2;
  color[2] = 7;
  color[3] = 6;
  color[4] = 5;
  color[5] = 5;
 // color[6] = 15;
 // color[7] = 1;
  int style[NCOMPMAX];
  style[0] = 1001;
  style[1] = 1001;
  style[2] = 1001;
  style[3] = 1001;
  style[4] = 1001;
  style[5] = 1001;
//  style[6] = 1001;
 // style[7] = 1001;
  double size[NCOMPMAX];
  int hitolo[NCOMPMAX];
  for (int i=0;i<nComponents;i++){
 //   hMC[isample][ibin][i][iatt]->SetLineColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillColor(color[i]);
    hMC[isample][ibin][i][iatt]->SetFillStyle(style[i]);
    size[i] = hMC[isample][ibin][i][iatt]->Integral();
    hitolo[i]=i;
  }
  int nswitch=1;
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
  hstack->Draw("h");
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
//  Leg->AddEntry(hMC[isample][ibin][5][iatt],"Single #pi^{+}","F");
 // Leg->AddEntry(hMC[isample][ibin][6][iatt],"Other","F");
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
  useSplineFlg=0;
  useNormFlg=1;
  return;
}


//constructor from a parameter file
histoManager::histoManager(const char* parfile, int nmode, bool separateneutmode)
  : nModes(nmode)
  , separateNeutMode(separateneutmode)
{
  
  //read in parameters
  sharedPars* runPars = new sharedPars(parfile);
  runPars->readParsFromFile();
 
  //set MC normalization
  normFactor = runPars->normFactor;

  //read in previously created histograms
  int nsamp = runPars->nSamples;
  int nbin  = runPars->nFVBins;
  int ncomp = runPars->nComponents;
  int natt  = runPars->nAttributes;
  readFromFile(runPars->hFactoryOutput.Data(),nsamp,nbin,ncomp,natt,nmode);

  //setup fit parameters
  fitPars = new atmFitPars(parfile);
  
  //get name string
  nameTag = runPars->globalRootName.Data();
 
  //setup histogram bin splines if using them 
  if (runPars->useSplinesFlg){
    readSplinesFromFile(runPars->splineFactoryOutput.Data());
  }
  
  // use norm pars by default
  useNormFlg =1.;
  return;
}

void histoManager::setFitPars(atmFitPars* thepars){

  // fix the pointer the the atmFitPars object that manages all of the necessary parameters
  fitPars=thepars;

  // use the normalization constant in the parameters
  normFactor = fitPars->normFactor;

  //
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

#ifdef T2K
void histoManager::setHistogram(int isamp, int ibin, int icomp, int imode, int iatt, int dataflg, TH1D *h)
{
  if (!dataflg){
    hMCNeut[isamp][ibin][icomp][imode][iatt] = h; 
  }
  else{
    hData[isamp][ibin][iatt] = h;
  }
}

void histoManager::setNominalHistogram(int isamp, int ibin, int icomp, int imode, int iatt, TH1D *h)
{
  hMCNeutNom[isamp][ibin][icomp][imode][iatt] = h;
}

void histoManager::fillHistogram(int isamp, int ibin, int icomp, int imode, int iatt, double value, double weight)
{
  hMCNeut[isamp][ibin][icomp][imode][iatt]->Fill(value,weight);
}

void histoManager::fillNominalHistogram(  int isamp, int ibin, int icomp, int imode, int iatt, double value, double weight)
{
  hMCNeutNom[isamp][ibin][icomp][imode][iatt]->Fill(value,weight);
}
#endif

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
void histoManager::readFromFile(const char* rootname,int nsamp,int nbin,int ncomp,int natt, int nmode){

  nSamples = nsamp;
  nBins    = nbin;
  nComponents = ncomp;
  nAttributes  = natt;
  nModes = nmode;
  TString filename = rootname;

  //file containing histograms 
  fin = new TFile(filename.Data());

  ///////////////////////////////////////////////
  //get normalization factor between data and MC
  //TH1D* htmp = (TH1D*)fin->Get("hnorm");

  ///////////////////////////////////////////////////
  //Set default sumw2 so errors are handled correctly
 // htmp->SetDefaultSumw2();

#ifdef T2K
  normFactor = SCALING;
#else
  //normFactor=htmp->GetBinContent(1);
#endif

  //////////////////////////////////////////////
  //read in histograms by name
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

  if (!separateNeutMode) {
    //setup mc histos
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
      	for (int icomp=0;icomp<nComponents;icomp++){
	        for (int iatt=0;iatt<nAttributes;iatt++){
	          hname = "hmc_";
	          hname.Append(Form("samp%d_bin%d_comp%d_att%d",isamp,ibin,icomp,iatt));
	          cout<<"Getting histogram: "<<hname.Data()<<endl;
	          hMC[isamp][ibin][icomp][iatt] = (TH1D*)fin->Get(hname.Data());
	          //set histogram mean array
//	          hMCMean[isamp][ibin][icomp][iatt] = hMC[isamp][ibin][icomp][iatt]->GetMean();
	          hMCMean[isamp][ibin][icomp][iatt] = hMC[isamp][ibin][icomp][iatt]->GetBinCenter(
                                                hMC[isamp][ibin][icomp][iatt]->GetMaximumBin());
	        }
      	}
      }
    }
  }

#ifdef T2K
  else {
    // setup mc histos
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
      	for (int icomp=0;icomp<nComponents;icomp++){
	        for (int imode = 0; imode < nModes; ++imode) {
	          for (int iatt=0;iatt<nAttributes;iatt++){
	            hname = "hmc_";
	            hname.Append(Form("samp%d_bin%d_comp%d_mode%d_att%d",isamp,ibin,icomp,imode,iatt));
	            //cout<<"Getting histogram: "<<hname.Data()<<endl;
	            hMCNeut[isamp][ibin][icomp][imode][iatt] = (TH1D*)fin->Get(hname.Data());
	            hname.Append("nom");
	            hMCNeutNom[isamp][ibin][icomp][imode][iatt] = (TH1D*)fin->Get(hname.Data());
	          }
	        }
	      }
      }
    }
  }
#endif

  ///////////////////////////////////////////////////////////////
  //initialize arrays of modified histograms using the histograms
  //from the file as templates
  initHistos(); 

  /////////////////////////////////////////////////////
  return;
}

///////////////////////////////////////////
//initializes arrays of modified histograms
void  histoManager::initHistos(){
  TH1::SetDefaultSumw2(kTRUE);
  TH1D* hsumset;
  hsumset->SetDefaultSumw2(kTRUE);
  cout<<"histoManager: initializing modified histograms..."<<endl;
  TString hmodname;
  if (!separateNeutMode) {
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
      	for (int iatt=0;iatt<nAttributes;iatt++){
      	  TString sumname = hMC[isamp][ibin][0][iatt]->GetName();
      	  sumname.Append("_sum");
      	  hSumHisto[isamp][ibin][iatt] = (TH1D*)hData[isamp][ibin][iatt]->Clone(sumname.Data());
	        sumname.Append("_modified");
      	  hSumHistoMod[isamp][ibin][iatt] = (TH1D*)hData[isamp][ibin][iatt]->Clone(sumname.Data()); 
          hSumHistoMod[isamp][ibin][iatt]->Reset();
	        for (int icomp=0;icomp<nComponents;icomp++){
	          hmodname = hMC[isamp][ibin][icomp][iatt]->GetName();
	          hmodname.Append("_mod");
	          hMCModified[isamp][ibin][icomp][iatt] = (TH1D*)hData[isamp][ibin][iatt]->Clone(hmodname.Data());
	          hMCModified[isamp][ibin][icomp][iatt]->Reset();
	       }
      	}
      }
    }
  }

#ifdef T2K
  else {
    for (int isamp=0;isamp<nSamples;isamp++){
      for (int ibin=0;ibin<nBins;ibin++){
      	for (int iatt=0;iatt<nAttributes;iatt++){
      	  TString sumname = hMCNeut[isamp][ibin][0][0][iatt]->GetName();
       	  sumname.Append("_sum");
       	  hSumHisto[isamp][ibin][iatt] = (TH1D*)hData[isamp][ibin][iatt]->Clone(sumname.Data());
	        sumname.Append("_modified");
	        hSumHistoMod[isamp][ibin][iatt] = (TH1D*)hData[isamp][ibin][iatt]->Clone(sumname.Data()); 
       	  for (int icomp=0;icomp<nComponents;icomp++){
	          for (int imode = 0; imode < nModes; ++imode) {
	            hmodname = hMCNeut[isamp][ibin][icomp][imode][iatt]->GetName();
	            hmodname.Append("_mod");
	            //std::cout<<hmodname<<std::endl;
	            hMCNeutModified[isamp][ibin][icomp][imode][iatt] = (TH1D*)hMC[isamp][ibin][icomp][iatt]->Clone(hmodname.Data());
	            hMCNeutModified[isamp][ibin][icomp][imode][iatt]->Reset();
	          }
	        }
	      }
      }
    }
  }
#endif
  return; 
}

void histoManager::readSplinesFromFile(const char* fname){

  //open spline parameter file
  TFile splineFile(fname);
  TTree* splinePars = (TTree*)splineFile.Get("splinePars");
  splineParReader* parReader = new splineParReader(splinePars);
  splinePars->GetEntry(0);
  int nsyspartot = parReader->nsyspartot;


  // make spline objects that will be filled with the spline parameters
  TString splinename;
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
  }
#ifdef T2K
  else {
    for (int ibin=0;ibin<nBins;ibin++){
      for (int isamp=0;isamp<nSamples;isamp++){
       	for (int icomp=0;icomp<nComponents;icomp++){
	        for (int imode = 0; imode < nModes; ++imode) {
	          for (int iatt=0;iatt<nAttributes;iatt++){
	            splinename="splinefor_";
	            splinename.Append(hMCNeut[isamp][ibin][icomp][imode][iatt]->GetName());
	            moreSplines[isamp][ibin][icomp][imode][iatt] =
                new hSplines(hMCNeutNom[isamp][ibin][icomp][imode][iatt],nsyspartot,splinename.Data());
	          }
	        }
	      }
      }
    }
  }
#endif

  //build the splines
  //double Y[parReader->npoints];
  //double X[parReader->npoints];
  cout<<"histoManager: Bulding histogram splines"<<endl;
  for (int ispline=0;ispline<splinePars->GetEntries();ispline++){
    splinePars->GetEntry(ispline);
    double Y[parReader->npoints];
    double X[parReader->npoints];
    for (int hbin=0;hbin<=parReader->nhistobins;hbin++){
      for (int ipt=0;ipt<parReader->npoints;ipt++){
        Y[ipt] = parReader->binWeight[ipt][hbin];
        X[ipt] = parReader->systParValues[ipt];
      }
      if (!separateNeutMode) {
      	theSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nattribute]
	       ->buildSpline(hbin,parReader->nsystpar,X,Y,parReader->npoints);
      }
#ifdef T2K
      else {
      	moreSplines[parReader->nsample][parReader->nbin][parReader->ncomponent][parReader->nmode][parReader->nattribute]
          ->buildSpline(hbin, parReader->nsystpar, X, Y, parReader->npoints);
      }
#endif
    }
  }
  splineFile.Close();
  useSplineFlg=1;
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this constructor creates a "blank" histogram manager that is used by histoFactory to 
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
  nComponents = 2;
  nAttributes = 1;
  nBins = 1;
  useSplineFlg=0; 
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
  //fitPars->norm = (double)nptsdata/(nptstot);
  normFactor = (double)nptsdata/(nptstot);
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

void histoManager::printBreakdownPlots(const char* plotdir){

  // directory where plots will be saved
  TString plotoutdir = plotdir;

  TString plotbasename = "mc_breakdown_";

  TCanvas* cc = new TCanvas("cc","cc",800,700);

  for (int ibin=0; ibin<nBins; ibin++){
    for (int isamp=0; isamp<nSamples; isamp++){
      for (int iatt=0; iatt<nAttributes; iatt++){
        TString tag = Form("comp_%d_bin_%d_att_%d.png",isamp,ibin,iatt);
        TString plotname = plotoutdir.Data();
        plotname.Append(plotbasename.Data());
        plotname.Append(tag.Data());
        showMCBreakdown(isamp,ibin,iatt);
        cc->Print(plotname.Data());
      }
    }
  }

  return;
}















#endif
