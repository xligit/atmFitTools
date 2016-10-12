{

 ///////////////////////////////////////
 //load class
 gROOT->ProcessLine(".L preProcess.cxx+");

 ///////////////////////////////////////
 //setup and run preprocessing object
 preProcess* preproc = new preProcess();
 preproc->setParFileName("atmparsE.dat");
 preproc->runPreProcessing();


}
