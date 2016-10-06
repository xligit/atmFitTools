{
gROOT->ProcessLine(".L calcResErr.C+");
gROOT->ProcessLine(".L histoTransforms.cxx+");

TH1D* h1 = testBumpD(100000,1,0,"h1");
h1->SetLineColor(kRed);
TH1D* h1c = (TH1D*)h1->Clone("h1c");
TH1D* h2 = testBumpD(105000,1.1,0.1,"h2");
TH1D* h2c = (TH1D*)h2->Clone("h2c");


TH1D* hb1 = testBumpD(100000,1,0,"hb1");
hb1->Rebin(2);
hb1->SetLineColor(kRed);
TH1D* hb1c = (TH1D*)h1->Clone("hb1c");
TH1D* hb2 = testBumpD(105000,1.1,0.1,"hb2");
hb2->Rebin(2);
TH1D* hb2c = (TH1D*)h2->Clone("hb2c");


TH1D* hd1 = testBumpD(1000,1,0,"d1");
hd1->SetLineColor(kRed);
TH1D* hd1c = (TH1D*)h1->Clone("hd1c");
TH1D* hd2 = testBumpD(1050,1.1,0.1,"d2");
TH1D* hd2c = (TH1D*)h2->Clone("hd2c");

h1->Draw();
h2->Draw("same");
}
