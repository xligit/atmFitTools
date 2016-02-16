{
 gROOT->ProcessLine(".L atmFitPars.C++"); 
// gROOT->ProcessLine(".L fQreader.C");
 gROOT->ProcessLine(".L postfitCalculator.C++");
// TFile f("nsk1.root");
// TTree* tt = (TTree*)f.Get("nsk");
 postfitCalculator *pcalc = new postfitCalculator("cosmicpars.dat");
// TChain chmc("h1");
// chmc.Add("./rootfiles/cosmic*ppmc*.root");
 //TChain chdat("h1");
 //chdat.Add("./rootfiles/cosmic*ppata*.root");
 //TTree* trmc = (TTree*)(&chmc);
 //TTree* trdata = (TTree*)(&chdat);
 //pcalc->setMCTree(trmc);
 //pcalc->NMCEvents=50000;
 //pcalc->cosmicPostFitAnalysis();
 //get previously made tree

// pcalc->setNSKTree(tt);
}
