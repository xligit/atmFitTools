{
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

  histoManager *hManager = new histoManager(hfilename.c_str(),pars->nSamples,pars->nBins,pars->nComponents,pars->nAttributes, pars->nModes, true);
  hManager->setFitPars(pars);
  hManager->readSplinesFromFile(splinefile.c_str(),pars->nSysPars);
  /*
    hSpline *h = hManager->moreSplines[0][2][1][0][0];
    h->draw2D(21, 53);
    h->evaluate();
   */
}
