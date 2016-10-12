
{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc= new histoCompare("atmparsE.dat");
// histoCompare* hc= new histoCompare("fakepars1.dat");

 // make fake data
 hc->setPar(1,50);
 hc->setPar(11,50);
 hc->setPar(21,50);
 hc->setPar(31,50);
 hc->setPar(41,50);
 hc->setPar(51,50);
 hc->makeFakeData();
// hc->setPar(1,0);
// hc->setPar(11,0);
// hc->setPar(21,0);
// hc->setPar(31,0);
// hc->setPar(41,0);
// hc->setPar(51,0);



// histoCompare* hc= new histoCompare("fakepars1.dat");
// histoCompare* hc= new histoCompare("hpi0pars.dat");
// hc->LnLFit();

}


