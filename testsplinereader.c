{
 TFile f("splineOut_debug.root");
 TTree* st = (TTree*)f.Get("splinePars");
 gROOT->ProcessLine(".L splineParReader.C+");
 splineParReader* reader = new splineParReader(st);
 st->GetEntry(0);


}
