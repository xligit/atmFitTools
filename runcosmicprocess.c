{

 ///////////////////////////////////////
 //load classes
 gROOT->ProcessLine(".L preProcess.C++");
// gROOT->ProcessLine(".L histoFactory.C++");
// gROOT->ProcessLine(".L splineFactory.C++");
 gROOT->ProcessLine(".x ~/style.c");

 preProcess* preproc = new preProcess();
 preproc->setParFileName("cosmicpars.dat");
 preproc->setWeightHistogram("./rootfiles/fine2.root","hpd");
 preproc->runPreProcessing();


}

