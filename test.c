{

 ///////////////////////////////////////
 //load class
 gROOT->ProcessLine(".L preProcess.cxx+");

 ///////////////////////////////////////
 //setup and run preprocessing object
 preProcess* preproc = new preProcess();
 preproc->setParFileName("atmparsE.dat");
 preproc->makeTestFiles("./rootfiles/",0,100000,10000,1);


}
