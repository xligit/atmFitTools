{

 ///////////////////////////////////////
 //load class
 gROOT->ProcessLine(".L preProcess.cxx+");

 ///////////////////////////////////////
 //setup and run preprocessing object
 preProcess* preproc = new preProcess();
 preproc->setParFileName("shimpars.dat");
 preproc->fakeShiftFlg = 0.;
 preproc->runPreProcessing();


}
