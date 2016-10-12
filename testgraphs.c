{

gROOT->ProcessLine(".L histoTransforms.cxx+");

TH1D* h = testBumpD(10000,5,5);
TH1D* hc = (TH1D*)h->Clone("clone");
hc->Scale(1./hc->GetBinWidth(1));
//TGraph* g = new TGraph(hc);
TGraph* g = histo2graph(hc);
shiftGraph(g,1.0,2.0);
double loss = graph2histo(g,hc);
double norm=h->Integral()-loss;
//hc->SetLineColor(kRed);
h->Draw();
hc->Scale(norm/hc->Integral());
hc->Draw("same");
}
