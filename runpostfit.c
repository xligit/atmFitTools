{
 gROOT->ProcessLine(".L atmFitPars.C++"); 
 gROOT->ProcessLine(".L fQreader.C");
 gROOT->ProcessLine(".L postfitCalculator.C++");
 postfitCalculator *pcalc = new postfitCalculator("mcmctree.root","sharedpars.dat");
 TChain chmc("h1");
 chmc.Add("./rootfiles/test1_MC*.root");
 TTree* trmc = (TTree*)(&chmc);
 pcalc->setMCTree(trmc);
 pcalc->nmcevents=10000;
}
