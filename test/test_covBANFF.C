void test_covBANFF()
{
  gROOT->Reset();
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C" );

  covBANFF *cov_prefit = new covBANFF("prefit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_joint_mecnubar_20151102.root", 1, false);
  //cov_prefit->PrintPars();
  cov_prefit->PrintNominal();

  covBANFF *cov_postfit = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_joint_mecnubar_20151102.root", 1, true);
  //cov_postfit->PrintPars();
  cov_postfit->PrintNominal();

}
