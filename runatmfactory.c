{

  ///////////////////////////////////////
  //load classes
  // gROOT->ProcessLine(".L histoManager.cxx++");
//  gROOT->ProcessLine(".L histoFactory.cxx+");
//  gROOT->ProcessLine(".L splineFactory.cxx+");
//  gROOT->ProcessLine(".L atmFitPars.cxx+"):
//  gROOT->ProcessLine(".x ~/style.c");

  //run histo factory from parameter file
  histoFactory* hfact = new histoFactory("atmpars1.dat");
  hfact->runHistoFactory();

  //run spline factory form parameter file
  //splineFactory *sfact = new splineFactory("atmpars1.dat"); //< atmospheric pars
  //sfact->runSplineFactory();


}
