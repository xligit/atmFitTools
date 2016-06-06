{
 gROOT->ProcessLine(".L sharedPars.cxx++");
 
 sharedPars* spar = new sharedPars("cosmicpars.dat");
 spar->readParsFromFile();

}
