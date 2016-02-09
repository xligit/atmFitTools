{

  ///////////////////////////////////////
  //load classes
  gROOT->ProcessLine(".L histoFactory.C++");
  gROOT->ProcessLine(".L splineFactory.C++");
  gROOT->ProcessLine(".x ~/style.c");

  //run histo factory from parameter file
  histoFactory* hfact = new histoFactory("cosmicpars.dat");
//  hfact->runHistoFactory();


  //run spline factory form parameter file
//  splineFactory *sfact = new splineFactory("sharedpars.dat"); //< atmospheric pars
//  splineFactory *sfact = new splineFactory("cosmicpars.dat"); //< cosmic pars
//  sfact->runSplineFactory();

/* //run histo factory

  //////////////////////////////////////////////////////
  //prefix for output files
  TString nametag = "nom2";
  TString factoryoutputname = "./rootfiles/";
  factoryoutputname.Append(nametag.Data());


  ///////////////////////////////////////////////////
  //set output file names
  TString hFactoryOutputName = "./rootfiles/histoFactory_fake1.root";
  TString sFactoryOutputName = "./rootfiles/splineFactory_fake1.root";


  ////////////////////////////////////////////////////
  //setup chains from post-processing files
  TString mcfiles= "./rootfiles/fake1_MC*.root";
  TString datfiles= "./rootfiles/fake1_Data_*.root";
  TChain* chmc=new TChain("h1");
  TChain *chdata=new TChain("h1");
  chmc->Add(mcfiles.Data());
  chdata->Add(datfiles.Data());
  TTree* trmc  = (TTree*)chmc;
  TTree* trdata = (TTree*)chdata;

  ///////////////////////////////////////////////
  //create and run histogram factory
  histoFactory* hfact = new histoFactory(3,3,7,factoryoutputname.Data()); 
  hfact->addAttribute(1);
  hfact->setDataTree(trdata);
  hfact->setMCTree(trmc);
  hfact->setOutputFileName(hFactoryOutputName.Data());
  hfact->init();
  hfact->fillHistos();
  hfact->normalizeHistos();
  hfact->saveToFile();
  TString factoryoutputname = hfact->getOutputFileName();
 
*/

   //spline factory
  ///////////////////////////////////////////////
/*  //create and run spline factory 
  TString splinename = "./rootfiles/";
  splinename.Append(nametag.Data());
  splineFactory* s = new splineFactory(3,3,7,1,1,splinename.Data());
  s->makeManagerFromFile(factoryoutputname.Data());
  s->setOutputFileName(sFactoryOutputName.Data());
  s->setupHistos();
  s->setupSystPars();
  s->setMCTree(trmc);
  s->buildTheSplines();
 */ 
}
