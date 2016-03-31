{
  //gROOT->ProcessLine(".L hSplines.C++");
  //gROOT->ProcessLine(".L histoManager.C++");
  //gROOT->ProcessLine(".L Tool_CompareToEventByEvent.C++");
  //gROOT->ProcessLine(".L atmFitPars.C++");
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");

  bool separateNeutMode = true;
  covXsec *cov = new covXsec("xsec_cov", "/home/xiaoyue/atmFitTools/rootfiles/xsec_covariance_2015a_q3_1.2_withRPA_v1.root");
  atmFitPars *pars = new atmFitPars("/home/xiaoyue/atmFitTools/atmpars.dat", cov);
  pars->printPars();

  double sigma = 1.0;
  int isys;
  for (int i = 0; i < pars->nSysPars; ++i) {
    if (pars->sysParName[i].find("pF_O")!=std::string::npos ||
	pars->sysParName[i].find("MAQE")!=std::string::npos ||
	pars->sysParName[i].find("MEC_O")!=std::string::npos ||
	pars->sysParName[i].find("EB_O")!=std::string::npos ||
	pars->sysParName[i].find("CA5")!=std::string::npos ||
	pars->sysParName[i].find("MANFFRES")!=std::string::npos ||
	pars->sysParName[i].find("BgRes")!=std::string::npos ||
	pars->sysParName[i].find("CCNUE_0")!=std::string::npos ||
	pars->sysParName[i].find("DISMPISHP")!=std::string::npos ||
	pars->sysParName[i].find("CCCOH_O_0")!=std::string::npos ||
	pars->sysParName[i].find("NCCOH_0")!=std::string::npos ||
	pars->sysParName[i].find("NCOTHER_0")!=std::string::npos ) {
      pars->setSysParameter(i, pars->sysParNom[i]+pars->sysParUnc[i]*sigma);
      //isys = i;
      //break;
      }
  }
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<", "<<pars->sysParName[i]<<" "<<pars->sysParNom[i]<<" "<<pars->sysPar[i]<<std::endl;
  }
  //pars->printPars();


  TChain* chmc = new TChain("h1");
  chmc->Add("/home/xiaoyue/atmFitTools/rootfiles/atmfit_ppmc_*_.root");

  compareToEventByEvent* comp = new compareToEventByEvent(pars,chmc,"/home/xiaoyue/atmFitTools/rootfiles/atmfit_histograms.root","/home/xiaoyue/atmFitTools/rootfiles/atmfit_splines.root", separateNeutMode);
  //comp->comparePrediction(0,0,0,0,0,(int)isys);//int isamp,int ibin, int icomp, int imode, int iatt, int isys);
  comp->comparePrediction(0,0,0,0,0,(bool)true);
  //std::cout<<pars->sysParName[isys]<<std::endl;
}
