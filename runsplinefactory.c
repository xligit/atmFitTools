{
 gROOT->ProcessLine(".L splineFactory.C+");
 TFile fmctree("nominal3.root");
 TTree* mctree = (TTree*)fmctree.Get("h1");
 splineFactory* s = new splineFactory(3,3,7,1,1,"debug");
 s->makeManagerFromFile("factoryOut_factorytest.root");
 //s->debugtest();
 s->setupHistos();
 s->setupSystPars();
 s->setMCTree(mctree);
 s->buildTheSplines();
}
