{
 gROOT->ProcessLine(".L atmFitPars.cxx++"); 
 gROOT->ProcessLine(".L postfitCalculator.cxx++");
 gROOT->ProcessLine(".x ~/style.c");
 postfitCalculator *pcalc = new postfitCalculator("pi0pars.dat");
 pcalc->attributeAnalysisFast();
 pcalc->calcPostFitHistos();

}
