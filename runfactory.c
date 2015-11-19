{
  gROOT->ProcessLine(".L histoFactory.C++");
  gROOT->ProcessLine(".L splineFactory.C++");
  gROOT->ProcessLine(".x ~/style.c");
  TString mcfiles= "./rootfiles/fake3_MC*.root";
  TString datfailes= "./rootfiles/nominal_Data*.root";
  TString nametag = "test1";

  TChain* chmc=new TChain("h1");
  TChain *chdata=new TChain("h1");

  chmc->Add(mcfiles.Data());
  chdata->Add(datfiles.Data());

  TTree* trmc  = (TTree*)chmc;
  TTree* trdata = (TTree*)chdata;

  TString factoryoutputname = "./rootfiles/";
  factoryoutputname.Append(nameTag.Data());
  histoFactory* hfact = new histoFactory(3,3,7,factoryoutputname.Data()); 
  hfact->addAttribute(1);
  hfact->setDataTree(trdata);
  hfact->setMCTree(trmc);
  hfact->init();
  //hm->readFromFile("hManager_atmos_emuratio");
  hfact->fillHistos();
  hfact->saveToFile();
  TString factoryoutputname = hfact->getOutputFileName();

  TString splinename="./rootfiles/";
  splinename.Append(namtag.Data());
  splineFactory* s = new splineFactory(3,3,7,1,1,splinename.Data());
  s->makeManagerFromFile(factoryoutputname.Data());
  s->setupHistos();
  s->setupSystPars();
  s->setMCTree(trmc);
  s->buildTheSplines();



}
