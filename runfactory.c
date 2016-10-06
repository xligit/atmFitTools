{

  ///////////////////////////////////////
  //load classes
 // gROOT->ProcessLine(".L histoManager.cxx++");
  gROOT->ProcessLine(".L histoFactory.cxx+");
  gROOT->ProcessLine(".L splineFactory.cxx+");
  gROOT->ProcessLine(".x ~/style.c");

////  run histo factory from parameter file
//  histoFactory* hfact = new histoFactory("fakepars3.dat");
  histoFactory* hfact = new histoFactory("shimpars.dat");
//  histoFactory* hfact = new histoFactory("fakepars2.dat");
//  histoFactory* hfact = new histoFactory("fakebumps.dat");
//  histoFactory* hfact = new histoFactory("atmparsE.dat");
//  histoFactory* hfact = new histoFactory("hpi0pars.dat");
  hfact->flgUseSample = 14228;
  hfact->runHistoFactory();
//



  //run spline factory form parameter file
 // splineFactory *sfact = new splineFactory("shimpars.dat"); //< atmospheric pars
//  sfact->buildSplineForPar(0);
 // sfact->runSplineFactory();


};
