{
 gROOT->ProcessLine(".L mcmcReader.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx++"); 
 gROOT->ProcessLine(".L postfitCalculator.cxx++");
 gROOT->ProcessLine(".x ~/style.c");
// postfitCalculator *pcalc = new postfitCalculator("atmparsE.dat");;
 postfitCalculator *pcalc = new postfitCalculator("hpi0pars.dat");;
 pcalc->attributeAnalysisFast();

 pcalc->calcPostFitHistos();

};
