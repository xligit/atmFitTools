{
 gROOT->ProcessLine(".L atmFitPars.C++"); 
 gROOT->ProcessLine(".L fQreader.C");
 gROOT->ProcessLine(".L postfitCalculator.C++");
// TFile f("nsk1.root");
// TTree* tt = (TTree*)f.Get("nsk");
 postfitCalculator *pcalc = new postfitCalculator("mcmctree.root","sharedpars.dat");
 TChain chmc("h1");
 chmc.Add("./rootfiles/test1_MC*.root");
 TTree* trmc = (TTree*)(&chmc);
 pcalc->setMCTree(trmc);
 pcalc->nmcevents=50000;

 //get previously made tree

// pcalc->setNSKTree(tt);
}
