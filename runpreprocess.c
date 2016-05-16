{

 ///////////////////////////////////////
 //load class
 gROOT->ProcessLine(".L preProcess.C+");

 ///////////////////////////////////////
 //setup and run preprocessing object
 preProcess* preproc = new preProcess();
 preproc->setParFileName("t2kmcpars.dat");
// preproc->runPreProcessing();


}
