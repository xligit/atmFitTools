{

  ///////////////////////////////////////
  //load classes
  // gROOT->ProcessLine(".L histoManager.C++");
//  gROOT->ProcessLine(".L histoFactory.C+");
//  gROOT->ProcessLine(".L splineFactory.C+");
//  gROOT->ProcessLine(".L atmFitPars.C+"):
//  gROOT->ProcessLine(".x ~/style.c");

  //run histo factory from parameter file
  histoFactory* hfact = new histoFactory("atmpars1.dat");
  hfact->runHistoFactory();

  //run spline factory form parameter file
  //splineFactory *sfact = new splineFactory("atmpars1.dat"); //< atmospheric pars
  //sfact->runSplineFactory();


}
