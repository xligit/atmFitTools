{
 gROOT->ProcessLine(".L histoManager.C++");
 gROOT->ProcessLine(".L hSplines.C++");
 gROOT->ProcessLine(".L histoCompare.C++");
 gROOT->ProcessLine(".L atmFitPars.C++");
 //gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc= new histoCompare("pi0pars.dat");
// hc->LnLFit();

}
