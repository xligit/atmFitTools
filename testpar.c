{
 gROOT->ProcessLine(".L sharedPars.C++");
 
 sharedPars* spar = new sharedPars("sharedpars.dat");
 spar->readParsFromFile();

}
