{
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");

  bool separateNeutMode = true;
  //covXsec *cov = new covXsec("xsec_cov", "/home/xiaoyue/atmFitTools/rootfiles/xsec_covariance_2015a_q3_1.2_withRPA_v1.root");
  covBANFF *cov = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools/rootfiles/postfit_data_1p1h_biascorrection_20160310.root", 1, true);
  atmFitPars *pars = new atmFitPars("/home/xiaoyue/atmFitTools/atmpars_mode.dat", cov);
  //pars->printPars();

  double sigma = 2;
  int isys = 53;

  TChain* chmc = new TChain("h1");
  chmc->Add("/home/xiaoyue/atmFitTools/rootfiles/neut_binning/atmfit_ppmc_*_.root");

  compareToEventByEvent* comp = new compareToEventByEvent(pars,chmc,"/home/xiaoyue/atmFitTools/rootfiles/neut_binning/atmfit_histograms.root","/home/xiaoyue/atmFitTools/rootfiles/neut_binning/atmfit_splines.root", separateNeutMode);

  for (int i = 0; i < pars->nSysPars; ++i) {
    if (pars->sysParName[i].find("pF_O")!=std::string::npos ||
	pars->sysParName[i].find("MAQE")!=std::string::npos ||
	pars->sysParName[i].find("MEC_O")!=std::string::npos ||
	pars->sysParName[i].find("EB_O")!=std::string::npos ||
	pars->sysParName[i].find("CA5")!=std::string::npos ||
	pars->sysParName[i].find("MANFFRES")!=std::string::npos ||
	pars->sysParName[i].find("BgRES")!=std::string::npos ||
	pars->sysParName[i].find("CCNUE")!=std::string::npos ||
	pars->sysParName[i].find("CCNUEBAR")!=std::string::npos ||
	pars->sysParName[i].find("MEC_NUBAR")!=std::string::npos ||
	pars->sysParName[i].find("DISMPISHP")!=std::string::npos ||
	pars->sysParName[i].find("CCCOH_O_0")!=std::string::npos ||
	pars->sysParName[i].find("NCCOH_0")!=std::string::npos ||
	pars->sysParName[i].find("NCOTHER_0")!=std::string::npos ||
	pars->sysParName[i].find("NC1GAMMA")!=std::string::npos ||
	pars->sysParName[i].find("FLUX_SUB")!=std::string::npos||
	pars->sysParName[i].find("FLUX_MUL")!=std::string::npos ) {
      pars->setSysParameter(i, pars->sysParNom[i]+pars->sysParUnc[i]*sigma);
    }
  }
  //for (int i = 0; i < pars->nSysPars; ++i) {
  //  std::cout<<i<<", "<<pars->sysParName[i]<<" "<<pars->sysParNom[i]<<" "<<pars->sysPar[i]<<std::endl;
  //}
  //comp->comparePrediction(0,2,0,0,0,(bool)true);

  for (int i = 0; i < pars->nSysPars; ++i) {
    if (pars->sysParName[i].find("MAQE")!=std::string::npos) {
      pars->setSysParameter(i, pars->sysParNom[i]+pars->sysParUnc[i]*sigma);
      isys = i;
    } else {
      pars->setSysParameter(i, pars->sysParNom[i]);
    }
  }
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<", "<<pars->sysParName[i]<<" "<<pars->sysParNom[i]<<" "<<pars->sysPar[i]<<std::endl;
  }
  std::cout<<pars->sysParName[isys]<<"\n"<<std::endl;
  for (int icom = 0; icom < pars->nComponents; ++icom) {
    comp->comparePrediction(0,2,icom,0,0,(int)isys);
    gROOT->ProcessLine(Form(".! mv compare53.pdf /home/xiaoyue/public_html/compare53_%d.pdf",icom));
  }

  return;

  for (int i = 0; i < pars->nSysPars; ++i) {
    if (pars->sysParName[i].find("CA5")!=std::string::npos) {
      pars->setSysParameter(i, pars->sysParNom[i]+pars->sysParUnc[i]*sigma);
      isys = i;
    } else {
      pars->setSysParameter(i, pars->sysParNom[i]);
    }
  }
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<", "<<pars->sysParName[i]<<" "<<pars->sysParNom[i]<<" "<<pars->sysPar[i]<<std::endl;
  }
  comp->comparePrediction(0,2,0,0,0,(int)isys);
  std::cout<<pars->sysParName[isys]<<"\n"<<std::endl;

  for (int i = 0; i < pars->nSysPars; ++i) {
    if (pars->sysParName[i].find("CCNUE_0")!=std::string::npos) {
      pars->setSysParameter(i, pars->sysParNom[i]+pars->sysParUnc[i]*sigma);
      isys = i;
    } else {
      pars->setSysParameter(i, pars->sysParNom[i]);
    }
  }
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<", "<<pars->sysParName[i]<<" "<<pars->sysParNom[i]<<" "<<pars->sysPar[i]<<std::endl;
  }
  comp->comparePrediction(0,2,0,0,0,(int)isys);
  std::cout<<pars->sysParName[isys]<<"\n"<<std::endl;

  for (int i = 0; i < pars->nSysPars; ++i) {
    if (pars->sysParName[i].find("FLUX_SUB")!=std::string::npos) {
      pars->setSysParameter(i, pars->sysParNom[i]+pars->sysParUnc[i]*sigma);
      isys = i;
    } else {
      pars->setSysParameter(i, pars->sysParNom[i]);
    }
  }
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<", "<<pars->sysParName[i]<<" "<<pars->sysParNom[i]<<" "<<pars->sysPar[i]<<std::endl;
  }
  comp->comparePrediction(0,2,0,0,0,(int)isys);
  std::cout<<pars->sysParName[isys]<<"\n"<<std::endl;

}
