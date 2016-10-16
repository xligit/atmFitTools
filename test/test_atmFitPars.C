void test_atmFitPars()
{
  gROOT->Reset();
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools_xligit/t2kLoadClasses.c");
  covBANFF *cov = new covBANFF("postfit_cov", "/home/xiaoyue/atmFitTools_xligit/rootfiles/postfit_data_1p1h_biascorrection_20160310.root", 1, true);
  /*
  std::cout<<"************** test covBANFF ******************"<<std::endl;
  for (int i=0; i< cov->getNPar(); ++i) {
    std::cout<<"# "<<i<<": "<<cov->getParName(i)<<", nominal "
	     <<cov->getNominal(i)<<", init "<<cov->getInit(i)<<", prior "
	     <<cov->getPrior(i)<<", proposed "<<cov->getProposed(i)
	     <<", current "<<cov->getCurrent(i)<<", uncertainty "
	     <<cov->getUncertainty(i)<<", lower "<<cov->getLow(i)<<", upper "
	     <<cov->getUp(i)<<std::endl;
  }
  std::cout<<"likelihood = "<<cov->getLikelihood()<<std::endl;
  cov->proposeStep();
  for (int i=0; i<cov->getNPar();++i) {
    std::cout<<"# "<<i<<": "<<cov->getParName(i)<<", initial "
	     <<cov->getInit(i)<<", current "<<cov->getCurrent(i)<<", proposed"
	     <<cov->getProposed(i)<<std::endl;
  }
  std::cout<<"likelihood = "<<cov->getLikelihood()<<std::endl;
  cov->acceptStep();
  for (int i=0; i<cov->getNPar();++i) {
    std::cout<<"# "<<i<<": "<<cov->getParName(i)<<", initial "
	     <<cov->getInit(i)<<", current "<<cov->getCurrent(i)<<", proposed"
	     <<cov->getProposed(i)<<std::endl;
  }
  std::cout<<"likelihood = "<<cov->getLikelihood()<<std::endl;
  for(int i=0; i<100; ++i) {
    cov->proposeStep();
    cov->acceptStep();
    std::cout<<"likelihood = "<<cov->getLikelihood()<<std::endl;
  }
  */
  atmFitPars *pars = new atmFitPars("/home/xiaoyue/atmFitTools_xligit/bin/atmpars.dat", cov);
  //pars->printParValues();
  //pars->printPars();
  std::cout<<"\n Initial value\n";
  /*
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<"th parameter "<<pars->sysParName[i]<<" has nominal value "<<pars->sysParNom[i]<<", and uncertain "<<pars->sysParUnc[i]<<", default "<<pars->parDefaultValue[pars->nTotPars - pars->nSysPars + i]
	     <<"; current "<<pars->sysPar[i]<<", proposed "<<pars->parsProp[pars->nTotPars - pars->nSysPars + i]<<std::endl;
  }
  std::cout<<"log prior "<<pars->calcLogPriors()<<std::endl;
  pars->proposeStep();
  std::cout<<"\nproposeStep"<<std::endl;
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<"th parameter "<<pars->sysParName[i]<<" has nominal value "<<pars->sysParNom[i]<<", and uncertain "<<pars->sysParUnc[i]<<", default "<<pars->parDefaultValue[pars->nTotPars - pars->nSysPars + i]
	     <<"; current "<<pars->sysPar[i]<<", proposed "<<pars->parsProp[pars->nTotPars - pars->nSysPars + i]<<std::endl;
  }
  std::cout<<"*** log "<<pars->calcLogPriors()<<std::endl;
  std::cout<<"\naccept step"<<std::endl;
  pars->acceptStep();
  for (int i = 0; i < pars->nSysPars; ++i) {
    std::cout<<i<<"th parameter "<<pars->sysParName[i]<<" has nominal value "<<pars->sysParNom[i]<<", and uncertain "<<pars->sysParUnc[i]<<", default "<<pars->parDefaultValue[pars->nTotPars - pars->nSysPars + i]
	     <<"; current "<<pars->sysPar[i]<<", proposed "<<pars->parsProp[pars->nTotPars - pars->nSysPars + i]<<std::endl;
  }
  std::cout<<"*** log\n "<<pars->calcLogPriors()<<std::endl;
  pars->proposeStep();
  pars->acceptStep();
  std::cout<<"*** log\n "<<pars->calcLogPriors()<<std::endl;
  */
  for (int i=0; i<pars->nTotPars; ++i) {
    std::cout<<i<<"th parameter "<<pars->parName[i]<<" has nominal value "<<pars->parDefaultValue[i]<<", and uncertain "<<pars->parUnc[i]<<", default "<<pars->parDefaultValue[i]
	     <<"; current "<<pars->pars[i]<<", proposed "<<pars->parsProp[i]<<std::endl;

  }
  std::cout<<"*** log\n "<<pars->calcLogPriors()<<std::endl;
  pars->proposeStep();
  std::cout<<"\n propose step"<<std::endl;
  for (int i=0; i<pars->nTotPars; ++i) {
    std::cout<<i<<"th parameter "<<pars->parName[i]<<" has nominal value "<<pars->parDefaultValue[i]<<", and uncertain "<<pars->parUnc[i]<<", default "<<pars->parDefaultValue[i]
	     <<"; current "<<pars->pars[i]<<", proposed "<<pars->parsProp[i]<<std::endl;

  }
  std::cout<<"*** log\n "<<pars->calcLogPriors()<<std::endl;
  pars->acceptStep();
  std::cout<<"\n accept step"<<std::endl;
  for (int i=0; i<pars->nTotPars; ++i) {
    std::cout<<i<<"th parameter "<<pars->parName[i]<<" has nominal value "<<pars->parDefaultValue[i]<<", and uncertain "<<pars->parUnc[i]<<", default "<<pars->parDefaultValue[i]
	     <<"; current "<<pars->pars[i]<<", proposed "<<pars->parsProp[i]<<std::endl;

  }
  std::cout<<"*** log\n"<<pars->calcLogPriors()<<std::endl;
  for(int i=0; i<10; ++i) {
    pars->proposeStep();
    pars->acceptStep();
    std::cout<<"*** log\n"<<pars->calcLogPriors()<<std::endl; 
  }
}
