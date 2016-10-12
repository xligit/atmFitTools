{

  ///////////////////////////////////////
  //load classes
 // gROOT->ProcessLine(".L histoManager.cxx++");
  gROOT->ProcessLine(".L histoFactory.cxx+");
  gROOT->ProcessLine(".L splineFactory.cxx+");
  gROOT->ProcessLine(".x ~/style.c");

////  run histo factory from parameter file
//  histoFactory* hfact = new histoFactory("fakepars1.dat");
//  histoFactory* hfact = new histoFactory("fakepars2.dat");
  histoFactory* hfact = new histoFactory("fakebumps.dat");
  hfact->runFakeFactory(9000,500,200,9000,500,200,50,1);
//  histoFactory* hfact = new histoFactory("atmparsE.dat");
//  histoFactory* hfact = new histoFactory("hpi0pars.dat");

//  hfact->runHistoFactory();
//



  //run spline factory form parameter file
//  splineFactory *sfact = new splineFactory("fakepars1.dat"); //< atmospheric pars
//  sfact->buildSplineForPar(0);
//  sfact->runSplineFactory();


};
