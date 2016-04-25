{
 gROOT->ProcessLine(".L atmFitPars.C++"); 
 gROOT->ProcessLine(".L postfitCalculator.C++");
 gROOT->ProcessLine(".x ~/style.c");
 postfitCalculator *pcalc = new postfitCalculator("pi0pars.dat");
 pcalc->attributeAnalysisFast();
 pcalc->calcPostFitHistos();

}
