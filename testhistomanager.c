{
gROOT->ProcessLine(".L histoManager.C++");
gROOT->ProcessLine(".x ~/style.c");
//TChain chdat("h1");
//TChain chmc("h1");
//chdat.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run073010.004*");
//chmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*001*fQ.root");
TFile fdata("postsift_atmospheric_FC_modmuoth.root");
TFile fmc("postsift_atmospheric_FC_2.root");
TTree* trdata = (TTree*)fdata.Get("h1");
TTree* trmc   = (TTree*)fmc.Get("h1");
histoManager* hm = new histoManager(3,1,7,"test1"); 
hm->addAttribute(1);
//hm->addAttribute(2);
hm->setDataTree(trdata);
hm->setMCTree(trmc);
hm->init();
//hm->readFromFile("hManager_atmos_emuratio");
hm->fillHistos();
hm->saveToFile();
}
