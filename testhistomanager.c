{
gROOT->ProcessLine(".L atmFitPars.cxx++");
gROOT->ProcessLine(".L hSplines.cxx++");
gROOT->ProcessLine(".L histoManager.cxx++");
//jgROOT->ProcessLine(".L hSplines.cxx++");



histoManager* hm = new histoManager("atmparsE.dat");

 

}
