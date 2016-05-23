{
  gROOT->ProcessLine(".x ~/rootlogon.C");
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");
  
  bool separateNeutMode = true;
  std::string card_file = "";
  std::string hfilename = "";
  std::string splinefile = "";
  if(separateNeutMode) {
    card_file = "/home/xiaoyue/atmFitTools/atmpars_mode.dat";
    hfilename = "/home/xiaoyue/atmFitTools/rootfiles/neut_binning/atmfit_histograms.root";
    splinefile = "/home/xiaoyue/atmFitTools/rootfiles/neut_binning/atmfit_splines.root";
  } else {
    card_file = "/home/xiaoyue/atmFitTools/atmpars.dat";
    hfilename = "/home/xiaoyue/atmFitTools/rootfiles/simple_binning/atmfit_histograms.root";
    splinefile = "/home/xiaoyue/atmFitTools/rootfiles/simple_binning/atmfit_splines.root";
  }

  covBANFF *cov = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_1p1h_biascorrection_20160310.root", 1, true);
  atmFitPars *pars = new atmFitPars(card_file.c_str(), cov);

  histoCompare* hc = new histoCompare(card_file.c_str(), separateNeutMode);
  hc->readFromFile(hfilename.c_str(), pars->nSamples, pars->nBins, pars->nComponents, pars->nModes, pars->nAttributes);
  hc->setupPars(pars); //setup parameters
  hc->setupSplines(splinefile.c_str(), pars->nSysPars);

// hc->setBinName(0,"bin0");
// hc->setBinName(1,"bin1");
// hc->setBinName(2,"bin2");
  hc->setCompName(0,"CC1e");
  hc->setCompName(1,"CC1mu");
  hc->setCompName(2,"CCeOth");
  hc->setCompName(3,"CCmuOth");
  hc->setCompName(4,"CCOth");
  hc->setCompName(5,"NCpi0");
  hc->setCompName(6,"NCOth");
  hc->setBinName(0,"FV0");
  hc->setBinName(1,"FV1");
  hc->setBinName(2,"FV2");
  hc->setBinName(3,"FV3");
  hc->setBinName(4,"FV4");
  hc->setBinName(5,"FV5");
  hc->setAttName(0,"emuPID");
  // hc->setRebinFactor(1);
  //hc->readFitPars("./rootfiles/fitpars_smooth_biasonly.root");
  //hc->flgFixAllSmearPars = 1;
  hc->tunePar=0.1;
  // hc->thePars->fixAllSmearPars(1);
  //hc->LnLFit();
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  hc->runMCMC(20000);
  sw->Stop();
  sw->Print("u");
  hc->saveFitPars("fitpars.root");
  // hc->addHistogram(hpidmc,0);
  // hc->addHistogram(hpidmcmu,0);
  // hc->addHistogram(hpiddat,1);
  // hc->addHistogram(hpidmcmu,0);
  // hc->addHistogram(hpiddatmu,1);
  
}
