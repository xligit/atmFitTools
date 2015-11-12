{
gROOT->ProcessLine(".L atmFitPars.C++");
gROOT->ProcessLine(".L hSplines.C++");
gROOT->ProcessLine(".L histoManager.C++");
//gROOT->ProcessLine(".x ~/style.c");
//gROOT->ProcessLine(".L atmFitPars.C++");

int nbin=3;
int ncomp=7;
int nsamp=3;
int natt=1;

//atmFitPars* fitpars = new atmFitpars(nbin,ncomp,nsamp,natt,0);
//TChain chdat("h1");
//TChain chmc("h1");
//chdat.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run073010.004*");
//chmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*001*fQ.root");
TChain *chdata = new TChain("h1");
TChain *chmc   = new TChain("h1");
//chdata->Add("./rootfiles/nominalRun_Data*.root");
//chmc->Add("./rootfiles/nominalRun_MC*.root");
chdata->Add("./rootfiles/multiSyst_Data*.root");
chmc->Add("./rootfiles/multiSyst_MC*.root");
atmFitPars* fitpars = new atmFitPars(nbin,ncomp,nsamp,natt,1);
 
TTree* trdata = (TTree*)chdata;
TTree* trmc   = (TTree*)trdata;

histoManager* hm = new histoManager("./rootfiles/multiSyst_hFactoryOutput.root",3,3,7,1); 
hm->readSplinesFromFile("./rootfiles/splineOutTest_splineOut.root",9);


hm->setFitPars(fitpars);
//hSplines* hs = hm->getSplines(0,0,0,0);

//hm->addAttribute(1);
//hm->addAttribute(2);
//hm->setDataTree(trdata);
//hm->setMCTree(trmc);
//hm->init();
//hm->readFromFile("hManager_atmos_emuratio");
//hm->fillHistos();
//hm->saveToFile();
}
