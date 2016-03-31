{
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");

  bool separateNeutMode = true;

  //run histo factory from parameter file
  //histoFactory* hfact = new histoFactory("/home/xiaoyue/atmFitTools/atmpars.dat", separateNeutMode);
  //hfact->runHistoFactory();

  covXsec *cov = new covXsec("xsec_cov", "/home/xiaoyue/atmFitTools/rootfiles/xsec_covariance_2015a_q3_1.2_withRPA_v1.root");
  atmFitPars *pars = new atmFitPars("/home/xiaoyue/atmFitTools/atmpars.dat", cov);
  //run spline factory form parameter file
  splineFactory *sfact = new splineFactory("/home/xiaoyue/atmFitTools/atmpars.dat", separateNeutMode); //< atmospheric pars
  sfact->setAtmFitPars(pars);
  /*
  sfact->makeManagerFromFile("/home/xiaoyue/atmFitTools/rootfiles/atmfit_histograms.root");
  sfact->setupHistos();
  sfact->setupSystPars();
  TChain c("h1");
  c.Add("/home/xiaoyue/atmFitTools/rootfiles/atmfit_ppmc_0_.root");
  sfact->setMCTree((TChain*)&c);
  sfact->mcTree->GetEntry(1);
  sfact->incrementSystPars(-2,3);
  std::cout<<sfact->mcEvt->nmode<<" "<<sfact->sysPar[3]<<" "<<sfact->mcEvt->evtweight<<std::endl;
  */
  sfact->runSplineFactory();
}
