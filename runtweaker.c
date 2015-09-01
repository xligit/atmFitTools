{
 gROOT->ProcessLine(".L tweaker.C++");
 TChain ch("h1");
 ch.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*fQ.root");
 tweaker* tw = new tweaker(&ch);
 tw->Alpha = 1.1;
 tw->Gamma = 5.2;
 tw->Beta  = 0.0;
 tw->meanDecayMom = 35.31;
 TH1D* hmom = new TH1D("hmom","hmom",50,0,80);
 tw->makeTweakTree("testy2");
}
