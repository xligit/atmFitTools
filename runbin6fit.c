
{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc= new histoCompare("fakepars1.dat");
// histoCompare* hc= new histoCompare("fakepars2.dat");

 // fix all other pars
 hc->thePars->fixAllPars(1);
 hc->thePars->fixPar[50] = 0;
 hc->thePars->fixPar[51] = 0;
//
// hc->thePars->fixPar[75] = 0;
 hc->thePars->fixPar[52] = 0;
 hc->thePars->fixPar[53] = 0;
 hc->thePars->fixPar[54] = 0;
 hc->thePars->fixPar[55] = 0;
 hc->thePars->fixPar[56] = 0;
 hc->thePars->fixPar[57] = 0;
 hc->thePars->fixPar[58] = 0;
 hc->thePars->fixPar[59] = 0;

 //hc->LnLFit();

}
