{
  gROOT->Reset();
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools_xligit/t2kLoadClasses.c");
  //  gDebug=2;  
  ///////////////////////////////////////
  //setup and run preprocessing object
  //preProcess* preproc1 = new preProcess();
  //preproc1->setParFileName("/home/xiaoyue/atmFitTools/atmpars.dat");
  //preproc1->runPreProcessing();
  
  t2kPreProcess* preproc2 = new t2kPreProcess();
  preproc2->setParFileName("/home/xiaoyue/atmFitTools_xligit/atmpar_mode.dat");
  preproc2->runPreProcessing();

}
