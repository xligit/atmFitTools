{
TFile f0("prof0.root");
TFile f1("prof1.root");

TH1F* h0 =(TH1F*)f0.Get("hprof");
TH1F* h1 =(TH1F*)f1.Get("hprof");

h0->SetLineColor(kRed);
h0->Draw();
h1->Draw("same");

}
