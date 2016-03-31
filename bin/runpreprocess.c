{
  gROOT->Reset();
  gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");
  
  ///////////////////////////////////////
  //setup and run preprocessing object
  preProcess* preproc = new preProcess();
  preproc->setParFileName("/home/xiaoyue/atmFitTools/atmpars.dat");
  preproc->runPreProcessing();
}
