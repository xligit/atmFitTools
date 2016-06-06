{

 ///////////////////////////////////////
 //load classes
 gROOT->ProcessLine(".L preProcess.cxx++");
// gROOT->ProcessLine(".L histoFactory.cxx++");
// gROOT->ProcessLine(".L splineFactory.cxx++");
 gROOT->ProcessLine(".x ~/style.c");

 preProcess* preproc = new preProcess();
 preproc->setParFileName("cosmicpars.dat");
 preproc->setWeightHistogram("./rootfiles/fine2.root","hpd");
 preproc->runPreProcessing();


}

