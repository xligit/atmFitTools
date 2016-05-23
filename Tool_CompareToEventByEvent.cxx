#include "Tool_CompareToEventByEvent.h"
#include "TCanvas.h"

void compareToEventByEvent::comparePrediction(int isamp,int ibin,int icomp,int iatt,int isys){
  //clone histogram
  TH1::SetDefaultSumw2(kTRUE);
  htmp1 = (TH1D*)hManager->getHistogram(isamp,ibin,icomp,iatt)->Clone("tmpclone");
  htmp1->Reset();
  htmp2 = hManager->getModHistogram(isamp,ibin,icomp,iatt);
  htmp3 = hManager->getHistogram(isamp,ibin,icomp,iatt);
  
  //fill histogram from tree
  for (int iev=0;iev<mcTree->GetEntries();iev++){
    mcTree->GetEvent(iev);
    fillAttributes();
    float weight = getEvtWeight(isys);
    std::cout<<"wgt = "<<std::endl;
    if ((mcEvt->nbin==ibin)&&(mcEvt->nsample==isamp)&&(mcEvt->ncomponent==icomp)){
      htmp1->Fill(att[iatt],weight);
    }
  }
  htmp1->SetLineColor(kBlack);
  htmp2->SetLineColor(kRed);
  htmp3->SetLineColor(kBlue);
  //htmp1->Draw("");
  //htmp1->Draw("samee");
  htmp2->Draw("samee");
  htmp2->Draw("sameh");
  htmp3->Draw("samee");
  htmp3->Draw("sameh");
  return;
}

void compareToEventByEvent::compareAllComp(int isamp, int ibin, int iatt, int isys)
{
  TH1::SetDefaultSumw2(kTRUE);
  htmp1 = (TH1D*)hManager->getHistogram(isamp,ibin,0,0,iatt)->Clone("tmpclone");
  htmp1->Reset();
  htmp2 = hManager->getSumHistogramMod(isamp,ibin,iatt,0);//getModHistogram(isamp,ibin,0,4,iatt);
  htmp3 = hManager->getSumHistogram(isamp,ibin,iatt,0);

  // fill histogram event-by-event

  for (int i = 0; i < mcTree->GetEntries(); ++i) {
    mcTree->GetEvent(i);
    fillAttributes();
    float weight = mcEvt->rfgweight * getEvtWeight(isys);
    if (i%50000==0) std::cout<<thePars->sysParName[isys]<<" "<<mcEvt->nmode<<" wgt = "<<weight<<std::endl;
    if ((mcEvt->nbin==ibin)&&(mcEvt->nsample==isamp)) {
      htmp1->Fill(att[iatt], weight);
    }
  }

  htmp1->SetLineColor(kBlack); // event-by-event modified histogram
  htmp1->SetLineWidth(2);
  htmp2->SetLineColor(kRed); // spline modified histogram
  htmp2->SetFillColor(kRed-9);
  htmp3->SetLineColor(kBlue); // original histogram
  TCanvas c;
  c.cd();
  htmp2->Draw("samee2");
  //htmp2->Draw("sameh");
  htmp1->Draw("sameh");
  htmp1->Draw("samee");
  htmp3->Draw("samee");
  htmp3->Draw("sameh");
  c.SaveAs(Form("compareAll%d.eps",isys));

}

void compareToEventByEvent::compareAllComp(int isamp, int ibin, int iatt)
{
  TH1::SetDefaultSumw2(kTRUE);
  htmp1 = (TH1D*)hManager->getHistogram(isamp,ibin,0,0,iatt)->Clone("tmpclone");
  htmp1->Reset();
  htmp2 = hManager->getSumHistogramMod(isamp,ibin,iatt,0);
  htmp3 = hManager->getSumHistogram(isamp,ibin,iatt,0);

  // fill histogram event-by-event
  for (int i = 0; i < mcTree->GetEntries(); ++i) {
    mcTree->GetEvent(i);
    fillAttributes();
    float weight = mcEvt->rfgweight;
    for (int isys = 0; isys < thePars->nSysPars; ++isys) {
      if (thePars->sysParName[isys].find("hc_")!=std::string::npos) continue;
      if (thePars->sysParName[isys].find("HAD")!=std::string::npos) continue;
      weight *= getEvtWeight(isys);
      if (i%50000==0) std::cout<<thePars->sysParName[isys]<<" "<<mcEvt->nmode<<" wgt = "<<weight<<std::endl;
    }
    if ((mcEvt->nbin==ibin)&&(mcEvt->nsample==isamp)) {
      htmp1->Fill(att[iatt], weight);
    }
  }

  htmp1->SetLineColor(kBlack); // event-by-event modified histogram
  htmp1->SetLineWidth(2);
  htmp2->SetLineColor(kRed); // spline modified histogram
  htmp2->SetFillColor(kRed-9);
  htmp3->SetLineColor(kBlue); // original histogram
  TCanvas c;
  c.cd();
  htmp2->Draw("samee2");
  //htmp2->Draw("sameh");
  htmp1->Draw("sameh");
  htmp1->Draw("samee");
  htmp3->Draw("samee");
  htmp3->Draw("sameh");
  c.SaveAs("compareAll.eps");

}

void compareToEventByEvent::comparePrediction(int isamp,int ibin,int icomp,int imode,int iatt,int isys){
  //clone histogram
  TH1::SetDefaultSumw2(kTRUE);
  htmp1 = (TH1D*)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->Clone("tmpclone");
  htmp1->Reset();
  htmp2 = hManager->getModHistogram(isamp,ibin,icomp,imode,iatt);
  htmp3 = hManager->getHistogram(isamp,ibin,icomp,imode,iatt);
  
  //fill histogram from tree
  for (int iev=0;iev<mcTree->GetEntries();iev++){
    mcTree->GetEvent(iev);
    fillAttributes();
    float weight = mcEvt->rfgweight * getEvtWeight(isys);
    if (iev%50000==0) std::cout<<thePars->sysParName[isys]<<" "<<mcEvt->nmode<<" wgt = "<<weight<<std::endl;
    if ((mcEvt->nbin==ibin)&&(mcEvt->nsample==isamp)&&(mcEvt->ncomponent==icomp)&&(mcEvt->nmode==imode)){
      htmp1->Fill(att[iatt],weight);
    }
  }
  htmp1->SetLineColor(kRed); // event-by-event modified histogram
  htmp1->SetLineWidth(2);
  htmp1->SetMarkerColor(kRed);
  htmp2->SetLineColor(kBlack); // spline modified histogram
  //htmp2->SetFillColor(kRed-9);
  htmp3->SetLineColor(kBlue); // original histogram
  std::cout<<"h1 "<<htmp1->Integral()<<"; h2 "<<htmp2->Integral()<<"; h3 "<<htmp3->Integral()<<std::endl;
  TCanvas c;
  c.cd();
  htmp1->Draw("h");
  htmp1->Draw("samee2");
  htmp2->Draw("samee");
  htmp2->Draw("sameh");
  htmp3->Draw("samee");
  htmp3->Draw("sameh");
  c.Print(Form("compare%d.pdf",isys));
  return;
}

void compareToEventByEvent::comparePrediction(int isamp,int ibin,int icomp,int imode,int iatt, bool all){
  //clone histogram
  TH1::SetDefaultSumw2(kTRUE);
  htmp1 = (TH1D*)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->Clone("tmpclone");
  htmp1->Reset();
  htmp2 = hManager->getModHistogram(isamp,ibin,icomp,imode,iatt);
  htmp3 = hManager->getHistogram(isamp,ibin,icomp,imode,iatt);
  
  //fill histogram from tree
  for (int iev=0;iev<mcTree->GetEntries();iev++){
    break;
    mcTree->GetEvent(iev);
    fillAttributes();
    float weight = mcEvt->rfgweight;
    for (int isys = 0; isys < thePars->nSysPars; ++isys) {
      weight *= getEvtWeight(isys);
      if (iev%50000==0) std::cout<<thePars->sysParName[isys]<<" "<<mcEvt->nmode<<" wgt = "<<weight<<std::endl;
    }
    if ((mcEvt->nbin==ibin)&&(mcEvt->nsample==isamp)&&(mcEvt->ncomponent==icomp)&&(mcEvt->nmode==imode)){
      htmp1->Fill(att[iatt],weight);
    }
  }
  htmp1->SetLineColor(kBlack);
  htmp1->SetLineWidth(2);
  htmp2->SetLineColor(kRed);
  htmp2->SetFillColor(kRed-9);
  htmp3->SetLineColor(kBlue);
  TCanvas c;
  c.cd();
  htmp2->Draw("samee2");
  //htmp2->Draw("sameh");
  htmp1->Draw("samee");
  htmp1->Draw("sameh");
  htmp3->Draw("samee");
  htmp3->Draw("sameh");
  c.SaveAs("compare.eps");
  return;
}

int compareToEventByEvent::getBest2RFitID(){
  int nfits = mcEvt->fqnmrfit;

  double ngLnLBest = 1000000.;
  int bestindex = 0;

  for (int ifit=0;ifit<nfits;ifit++){
    int fitID = mcEvt->fqmrifit[ifit]; //< fit fit ID code
    if ((fitID-320000000)<0) continue; //< we want best 2R fits
    if (mcEvt->fqmrnll[ifit]<ngLnLBest) bestindex = ifit;
  }
  return bestindex;
}

void compareToEventByEvent::fillAttributes(){
  //att[0] = mcEvt->fq1rnll[0][2]-mcEvt->fq1rnll[0][1];
  //att[1] = mcEvt->fq1rnll[1][2]-mcEvt->fq1rnll[1][1];
  //int ibest = getBest2RFitID();
  //double best1Rnglnl = fmin(mcEvt->fq1rnll[0][1],mcEvt->fq1rnll[0][2]);
  //att[1] = best1Rnglnl-mcEvt->fqmrnll[ibest];
  att[0] = mcEvt->attribute[0];
  att[1] = mcEvt->attribute[1];
}

float compareToEventByEvent::getEvtWeight(int ipar){
//  float ww = 1.;
  float ww = 1;//mcEvt->evtweight;
  int absmode = TMath::Abs(mcEvt->mode);
  float Enu     = mcEvt->pmomv[0];
  int  nutype  = mcEvt->ipnu[0];
  double evis = mcEvt->fq1rmom[0][1]; 
  int nmode = mcEvt->nmode; 
  if (thePars->sysParName[ipar].find("_C")!=std::string::npos ||
      thePars->sysParName[ipar].find("hc_")!=std::string::npos) {
    return ww;
  }
  if (thePars->sysParName[ipar].find("MAQE")!=std::string::npos) {
    if (nmode == 0) {
      ww *= mcEvt->byEv_maqe_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<mcEvt->byEv_maqe_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("pF_O")!=std::string::npos) {
    if (nmode == 0) {
      ww *= mcEvt->byEv_pfo_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<mcEvt->byEv_pfo_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("EB_O")!=std::string::npos) {
    if (nmode == 0) {
      ww *= mcEvt->byEv_ebo_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S"); 
      //std::cout<<thePars->sysParName[ipar]<<" "<<mcEvt->byEv_ebo_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("CA5")!=std::string::npos) {
    if (nmode == 1) {
      ww *= mcEvt->byEv_ca5_cc1pi_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      ////std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_ca5_cc1pi_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
    else if (nmode == 4) {
      ww *= mcEvt->byEv_ca5_ncpiz_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_ca5_ncpiz_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
    else if (nmode == 5) {
      ww *= mcEvt->byEv_ca5_ncpipm_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_ca5_ncpipm_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("MANFF")!=std::string::npos) {
    if (nmode == 1) {
      ww *= mcEvt->byEv_manff_cc1pi_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_manff_cc1pi_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
    else if (nmode == 4) {
      ww *= mcEvt->byEv_manff_ncpiz_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_manff_ncpiz_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
    else if (nmode == 5) {
      ww *= mcEvt->byEv_manff_ncpipm_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_manff_ncpipm_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("BgRES")!=std::string::npos) {
    if (nmode == 1) {
      ww *= mcEvt->byEv_bgscl_cc1pi_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_bgscl_cc1pi_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
    else if (nmode == 4) {
      ww *= mcEvt->byEv_bgscl_ncpiz_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_bgscl_ncpiz_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
    else if (nmode == 5) {
      ww *= mcEvt->byEv_bgscl_ncpipm_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_bgscl_ncpipm_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("DISMPISHP")!=std::string::npos) { 
    if (nmode == 3) {
      ww *= mcEvt->byEv_dismpishp_ccoth_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S");
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_dismpishp_ccoth_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]), 0, "S")<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("FLUX_SUB")!=std::string::npos) {
    if (evis<1300) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("FLUX_MUL")!=std::string::npos) { 
    if (evis>1300) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("HAD_MUL")!=std::string::npos) {}//do nothing for now
  /*else if (thePars->sysParName[ipar].find("RPA_O")!=std::string::npos) {
    if (nmode == 0) {
      ww *= mcEvt->byEv_rpa_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]));
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_rpa_ccqe_gr->Eval((thePars->sysPar[ipar]-thePars->sysParNom[ipar]))<<" | ";
    }
    }*/
  else if (thePars->sysParName[ipar].find("MEC_O")!=std::string::npos) {
    if (nmode == 8) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("MEC_NUBAR")!=std::string::npos) {
    if (nmode == 8 && nutype < 0) ww *= thePars->sysPar[ipar];
  }
  else if (thePars->sysParName[ipar].find("CCNUE")!=std::string::npos) {
    if (nutype == 12 && nmode <4) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("CCNUEBAR")!=std::string::npos) {
    if (nutype==-12 && nmode < 4) ww *= thePars->sysPar[ipar];
  }
  else if (thePars->sysParName[ipar].find("CCCOH_O_0")!=std::string::npos) {
    if (nmode == 3) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("NCOTHER_0")!=std::string::npos) {
    if (nmode == 7) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("NCCOH_0")!=std::string::npos) {
    if (nmode == 6) {
      ww *= thePars->sysPar[ipar];
      //std::cout<<thePars->sysParName[ipar]<<" "<<nmode<<" "<<thePars->sysPar[ipar]<<" | ";
    }
  }
  else if (thePars->sysParName[ipar].find("NC1GAMMA_O")!=std::string::npos) {
    if (nmode == 9) ww *= thePars->sysPar[ipar];
  }
  else if (thePars->sysParName[ipar].find("FSI_INEL_LO_E")!=std::string::npos) {}//do nothing for now
  else if (thePars->sysParName[ipar].find("FSI_INEL_HI_E")!=std::string::npos) {}//do nothing for now
  else if (thePars->sysParName[ipar].find("FSI_PI_PROD")!=std::string::npos) {}//do nothing for now
  else if (thePars->sysParName[ipar].find("FSI_PI_ABS")!=std::string::npos) {}//do nothing for now
  else if (thePars->sysParName[ipar].find("FSI_CEX_LO_E")!=std::string::npos) {}//do nothing for now
  else if (thePars->sysParName[ipar].find("FSI_CEX_HI_E")!=std::string::npos) {}//do nothing for now
  else {} // all other systematic parameters don't change the weights
  
  //std::cout<<"----------- "<<"mode="<<nmode<<" "<<ww<<" -----------"<<std::endl;
  return ww>0 ? ww:0;
}

compareToEventByEvent::compareToEventByEvent(atmFitPars* thepars, TChain* tr,const char* hfilename, const char* sfilename, bool separateNeutMode)
  : separateNeutMode(separateNeutMode)
{
  //setup mc tree
  mcTree = tr;
  mcEvt = new fQreader(mcTree);
  //setup fit pars
  thePars = thepars;
  //creaate manager
  if (!separateNeutMode) {
    hManager = new histoManager(hfilename,thePars->nSamples,thePars->nBins,thePars->nComponents,thePars->nAttributes);
  } else {
    hManager = new histoManager(hfilename,thePars->nSamples,thePars->nBins,thePars->nComponents,thePars->nAttributes, thePars->nModes, true);
  }
  //read in spline info
  hManager->setFitPars(thePars);
  hManager->readSplinesFromFile(sfilename,thePars->nSysPars); 
}


