{
  gROOT->ProcessLine(".L histoFactory.C++");
  gROOT->ProcessLine(".L splineFactory.C++");
  gROOT->ProcessLine(".x ~/style.c");
  TString mcfiles= "./rootfiles/nominal2_MC*.root";
  TString datfiles= "./rootfiles/nominal2_Data_*.root";
  TString nametag = "nom2";

  TChain* chmc=new TChain("h1");
  TChain *chdata=new TChain("h1");

  chmc->Add(mcfiles.Data());
  chdata->Add(datfiles.Data());

  TTree* trmc  = (TTree*)chmc;
  TTree* trdata = (TTree*)chdata;

  TString factoryoutputname = "./rootfiles/";
  factoryoutputname.Append(nametag.Data());
  histoFactory* hfact = new histoFactory(3,3,7,factoryoutputname.Data()); 
  hfact->addAttribute(1);
  hfact->setDataTree(trdata);
  hfact->setMCTree(trmc);
  hfact->setOutputFileName("./rootfiles/nom2_factoryOutput.root");
  hfact->init();
  //hm->readFromFile("hManager_atmos_emuratio");
  hfact->fillHistos();
  hfact->normalizeHistos();
  hfact->saveToFile();
  TString factoryoutputname = hfact->getOutputFileName();
 
  TString splinename = "./rootfiles/";
  splinename.Append(nametag.Data());
  splineFactory* s = new splineFactory(3,3,7,1,1,splinename.Data());
  s->makeManagerFromFile(factoryoutputname.Data());
  s->setOutputFileName("./rootfiles/nom2_splineFactoryOut.root");
  s->setupHistos();
  s->setupSystPars();
  s->setMCTree(trmc);
  s->buildTheSplines();

}
