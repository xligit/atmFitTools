{
 gROOT->ProcessLine(".L atmFitPars.C++"); 
 gROOT->ProcessLine(".L postfitCalculator.C++");
 postfitCalculator *pcalc = new postfitCalculator("mcmcsmooth.root","sharedpars.dat");
 TChain chmc("h1");
 chmc.Add("./rootfiles/nominal3_MC*.root");
 TTree* trmc = (TTree*)(&chmc);
 pcalc->setMCTree(trmc);
}
