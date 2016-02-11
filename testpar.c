{
 gROOT->ProcessLine(".L sharedPars.C++");
 
 sharedPars* spar = new sharedPars("cosmicpars.dat");
 spar->readParsFromFile();

}
