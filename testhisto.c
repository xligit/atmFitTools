{
 gROOT->ProcessLine(".L TH2FV.cxx");
 TH2FV* h = new TH2FV("h",4);
 h->Draw();
}
