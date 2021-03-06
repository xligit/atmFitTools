{
 gROOT->ProcessLine(".x ~/rootlogon.C");
 gROOT->ProcessLine(".x /home/xiaoyue/atmFitTools/loadClasses.C");
 TH1F* hpdat = new TH1F("hpdat","hpdat",50,0,100);
 TH1F* hpmc = new TH1F("hpmc","hpmc",50,0,100);
 TH1F* hpiddat = new TH1F("hpiddat","hpiddat",50,-1000,600);
 TH1F* hpidmc = new TH1F("hpidmc","hpidmc",50,-1000,600);
 TH1F* hpiddatmu = new TH1F("hpiddatmu","hpiddatmu",50,-1000,600);
 TH1F* hpidmcmu = new TH1F("hpidmcmu","hpidmcmu",50,-1000,600);

 histoCompare* hc = new histoCompare("comptest");
// hc->readFromFile("histos_test2",3,3,7,1);  //reads in histogram manager
 hc->readFromFile("./rootfiles/histoFactory_fake1.root",3,3,7,1);
// hc->readFromFile(".root",3,3,7,1);

 hc->setupPars(1); //setup parameters
 hc->thePars->setSysParUnc(0,0.05);
 hc->setupSplines("./rootfiles/splineFactory_fake1.root",9);
// hc->setupSplines("./rootfiles/splineOutTest_splineOut.root");
// hc->setupSplines("./rootfiles/nominalRun_splineOut.root");

// hc->hManager->useSplineFlg=1;
// hc->getTotSumSq();
// hc->setBinName(0,"bin0");
// hc->setBinName(1,"bin1");
// hc->setBinName(2,"bin2");
 hc->setCompName(0,"CC1e");
 hc->setCompName(1,"CC1mu");
 hc->setCompName(2,"CCeOth");
 hc->setCompName(3,"CCmuOth");
 hc->setCompName(4,"CCOth");
 hc->setCompName(5,"NCpi0");
 hc->setCompName(6,"NCOth");
 hc->setBinName(0,"FV0");
 hc->setBinName(1,"FV1");
 hc->setBinName(2,"FV2");
 hc->setAttName(0,"emuPID");
 hc->setAttName(0,"Other");
// hc->setRebinFactor(1);
 //hc->readFitPars("./rootfiles/fitpars_smooth_biasonly.root");
 //hc->flgFixAllSmearPars = 1;
 hc->tunePar=0.1;
// hc->thePars->fixAllSmearPars(1);
 //hc->LnLFit();
 hc->runMCMC(1000);
 //hc->saveFitPars("./rootfiles/fitpars.root");
// hc->addHistogram(hpidmc,0);
// hc->addHistogram(hpidmcmu,0);
// hc->addHistogram(hpiddat,1);
// hc->addHistogram(hpidmcmu,0);
// hc->addHistogram(hpiddatmu,1);

}
