void runfactory(bool i=true)
{
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");

  bool separateNeutMode = i;
  std::string card_file = "";
  if(separateNeutMode) card_file = "/home/xiaoyue/atmFitTools/atmpars_mode.dat";
  else card_file = "/home/xiaoyue/atmFitTools/atmpars.dat";

  //run histo factory from parameter file
  //histoFactory* hfact = new histoFactory(card_file.c_str(), separateNeutMode);
  //hfact->runHistoFactory();

  //covXsec *cov = new covXsec("xsec_cov", "/home/xiaoyue/atmFitTools/rootfiles/xsec_covariance_2015a_q3_1.2_withRPA_v1.root");
  covBANFF *cov = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_1p1h_biascorrection_20160310.root", 1, true);
  atmFitPars *pars = new atmFitPars(card_file.c_str(), cov);
  //run spline factory form parameter file
  splineFactory *sfact = new splineFactory(card_file.c_str(), separateNeutMode); //< atmospheric pars
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
