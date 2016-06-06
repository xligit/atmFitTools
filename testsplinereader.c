{
 TFile f("./rootfiles/multiSyst_splineOut.root");
 TTree* st = (TTree*)f.Get("splinePars");
 gROOT->ProcessLine(".L splineParReader.cxx+");
 splineParReader* reader = new splineParReader(st);
 st->GetEntry(0);
}
