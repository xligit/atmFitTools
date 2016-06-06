{
 gROOT->ProcessLine(".L histoManager.cxx++");
 gROOT->ProcessLine(".L hSplines.cxx++");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx++");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc= new histoCompare("cosmicpars.dat");
// hc->LnLFit();

}

