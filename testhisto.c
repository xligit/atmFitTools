{
 gROOT->ProcessLine(".L TH2FV.C");
 TH2FV* h = new TH2FV("h",4);
 h->Draw();
}
