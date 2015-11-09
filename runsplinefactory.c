{
 gROOT->ProcessLine(".L splineFactory.C+");

 TChain* chmc = new TChain("h1");
 chmc->Add("./rootfiles/nominalRun_MC*.root");
 TTree* trmc = (TTree*)chmc;
 TString splinename="./rootfiles/";
 splinename.Append("splineOutTest");
 splineFactory* s = new splineFactory(3,3,7,1,1,splinename.Data());
 s->makeManagerFromFile("./rootfiles/nominalRun_hFactoryOutput.root");
 s->setupHistos();
 s->setupSystPars();
 s->setMCTree(trmc);
 s->buildTheSplines();


 /*TFile fmctree("nominal3.root");
 TTree* mctree = (TTree*)fmctree.Get("h1");
 splineFactory* s = new splineFactory(3,3,7,1,1,"debug");
 s->makeManagerFromFile("factoryOut_factorytest.root");
 //s->debugtest();
 s->setupHistos();
 s->setupSystPars();
 s->setMCTree(mctree);
 s->buildTheSplines();
*/
}
