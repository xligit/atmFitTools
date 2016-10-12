
{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc= new histoCompare("fakepars1.dat");

 // setup priors
 double pi0bias = 10.;
 double pi0smear = 0.05;

 hc->thePars->setGausPrior(18,pi0smear);
 hc->thePars->setGausPrior(19,pi0bias);

 hc->thePars->setGausPrior(28,pi0smear);
 hc->thePars->setGausPrior(29,pi0bias);

 hc->thePars->setGausPrior(38,pi0smear);
 hc->thePars->setGausPrior(39,pi0bias);

 hc->thePars->setGausPrior(48,pi0smear);
 hc->thePars->setGausPrior(49,pi0bias);

 hc->thePars->setGausPrior(58,pi0smear);
 hc->thePars->setGausPrior(59,pi0bias);

 hc->thePars->setGausPrior(0,pi0smear);
 hc->thePars->setGausPrior(1,pi0bias);

}


