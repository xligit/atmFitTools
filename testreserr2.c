{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L calcResErr.C+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");

// histoCompare* hc= new histoCompare("atmparsE.dat");
// histoCompare* hc= new histoCompare("fakepars1.dat");
 histoCompare* hc= new histoCompare("hpi0pars.dat");
 hc->showFitResult(0,0,0);

// hc->LnLFit();

}


