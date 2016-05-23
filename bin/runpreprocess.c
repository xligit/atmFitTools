{
  gROOT->Reset();
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");
  
  ///////////////////////////////////////
  //setup and run preprocessing object
  //preProcess* preproc1 = new preProcess();
  //preproc1->setParFileName("/home/xiaoyue/atmFitTools/atmpars.dat");
  //preproc1->runPreProcessing();
  
  preProcess* preproc2 = new preProcess();
  preproc2->setParFileName("/home/xiaoyue/atmFitTools/atmpars_mode.dat");
  preproc2->runPreProcessing();

}
