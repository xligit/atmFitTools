{

 ///////////////////////////////////////
 //load class
 gROOT->ProcessLine(".L preProcess.cxx+");

 ///////////////////////////////////////
 //setup and run preprocessing object
 preProcess* preproc = new preProcess();
 preproc->setParFileName("fakepars3.dat");
// preproc->setParFileName("hpi0pars.dat");
// preproc->setParFileName("atmparsE.dat");
 preproc->fakeShiftFlg = 0.;
// preproc->fakeNormFlg = 0.;
 preproc->runPreProcessing();


}
