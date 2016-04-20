{
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");

  bool separateNeutMode = true;

  covBANFF *cov = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_joint_mecnubar_20151102.root", 1, true);
  atmFitPars *pars = new atmFitPars("/home/xiaoyue/atmFitTools/atmpars.dat", cov);

  const char * hfilename = "/home/xiaoyue/atmFitTools/rootfiles/atmfit_histograms.root";

  histoManager *hManager = new histoManager(hfilename,pars->nSamples,pars->nBins,pars->nComponents,pars->nAttributes, pars->nModes, true);
  hManager->setFitPars(pars);
  hManager->readSplinesFromFile("/home/xiaoyue/atmFitTools/rootfiles/atmfit_splines.root",pars->nSysPars);

}
