{
gROOT->ProcessLine(".L histoFactory.C++");
gROOT->ProcessLine(".x ~/style.c");
//TChain chdat("h1");
//TChain chmc("h1");
//chdat.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run073010.004*");
//chmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*001*fQ.root");
TFile fdata("nominal4.root");
//TFile fmc("nominal2.root");
TTree* trdata = (TTree*)fdata.Get("h1");
TChain* chmc = new TChain("h1");
chmc->Add("test_*.root");
TTree* trmc = (TTree*)(chmc);
//TTree* trmc   = (TTree*)fmc.Get("h1");
histoFactory* hfact = new histoFactory(3,3,7,"factorytest2"); 
hfact->addAttribute(1);
hfact->setDataTree(trdata);
hfact->setMCTree(trmc);
hfact->init();
//hm->readFromFile("hManager_atmos_emuratio");
hfact->fillHistos();
hfact->saveToFile();
}
