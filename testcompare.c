{
 gROOT->ProcessLine(".L histoManager.cxx++");
 gROOT->ProcessLine(".L hSplines.cxx++");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx++");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc = new histoCompare("comptest");
 histoManager* hm = new histoManager(800,700);
 hc->initialize(hm,hm->fitPars);
 hc->profileL(1,3,100,0);
// hc->useLnLType = 0;
// hc->cScale = 1.;
// hc->profileL(1,3,100,1);


}
