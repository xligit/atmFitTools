{
 gROOT->ProcessLine(".L splineFactory.C+");
 splineFactory* s = new splineFactory(3,3,7,1,1,"debug");
 s->makeManagerFromFile("factoryOut_factorytest.root");
 //s->debugtest();
 s->setupHistos();
}
