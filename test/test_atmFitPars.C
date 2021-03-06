void test_atmFitPars()
{
  gROOT->Reset();
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C" );
  //covXsec *cov = new covXsec("xsec_cov", "/home/xiaoyue/atmFitTools/rootfiles/xsec_covariance_2015a_q3_1.2_withRPA_v1.root");
  //cov->PrintPars();
  //cov->PrintNominal();
  covBANFF *cov = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_1p1h_biascorrection_20160310.root", 1, true);
  //cov_postfit->PrintPars();
  cov->PrintNominal();
  atmFitPars *pars = new atmFitPars("/home/xiaoyue/atmFitTools/atmpars_mode.dat", cov);
  //pars->printParValues();
  //pars->printPars();
  for (int i = 0; i < pars->nSysPars; ++i) {
    //std::cout<<i<<"th parameter "<<pars->sysParName[i]<<" has nominal value "<<pars->sysParNom[i]<<", and uncertain "<<pars->sysParUnc[i]<<std::endl;
  }

  TChain *c = new TChain("h1");
  c->Add("/home/xiaoyue/atmFitTools/rootfiles/neut_binning/atmfit_ppmc_0_.root");
  fQreader *r = new fQreader(c);
  r->FillMap();

  std::cout<<"\n------------------------------ propose step ----------------------"<<std::endl;
  pars->proposeStep();
  pars->printPars();
  std::cout<<"\n------------------------------ accept step -----------------------"<<std::endl;
  pars->acceptStep();
  pars->printPars();
}
