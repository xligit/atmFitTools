{

 ///////////////////////////////////////
 //load classes
 gROOT->ProcessLine(".L preProcess.C++");
 gROOT->ProcessLine(".x ~/style.c");

 preProcess* preproc = new preProcess();
 preproc->setParFileName("atmpars1.dat");
 // preproc->setWeightHistogram("./rootfiles/fine2.root","hpd");
 preproc->runPreProcessing();


}

