{



TLine* li = new TLine(50,0,50,100000);
li->SetLineWidth(3);
li->SetLineColor(kOrange);


// mcmc chains
TFile* ffit = new TFile("./mcmc/fake1_biasonly_minfit.root");
TFile* ftru = new TFile("./mcmc/fake1_biasonly_truefit.root");

TTree* tfit = (TTree*)ffit->Get("MCMCpath");;
TTree* ttru = (TTree*)ftru->Get("MCMCpath");;

double norm = (double)ttru->GetEntries()/(double)tfit->GetEntries();

// bin 0
//TH1D* h2 = new TH1D("hf","hf",100,0,90);
//TH1D* h1 = new TH1D("ht","ht",100,0,90);
tfit->Draw("par[25]>>hf");
ttru->Draw("par[25]>>ht","step>300000");
ht->Scale(hf->Integral()/ht->Integral());
ht->SetLineColor(kBlue);
ht->Draw();
hf->Draw("same");
li->Draw("same");

/*
// bin 1
tfit->Draw("par[25]>>hf");
ttru->Draw("par[25]>>ht","step>300000");
ht->Scale(hf->Integral()/ht->Integral());
ht->SetLineColor(kBlue);
ht->Draw();
hf->Draw("same");
li->Draw("same");
c1->Print("~/transfer/stuck_compare_bin0.png");
h2->Reset();
h1->Reset();

// bin 2
tfit->Draw("par[25]>>hf");
ttru->Draw("par[25]>>ht","step>300000");
ht->Scale(hf->Integral()/ht->Integral());
ht->SetLineColor(kBlue);
ht->Draw();
hf->Draw("same");
li->Draw("same");
h2->Reset();
h1->Reset();

// bin 3
tfit->Draw("par[25]>>hf");
ttru->Draw("par[25]>>ht","step>300000");
ht->Scale(hf->Integral()/ht->Integral());
ht->SetLineColor(kBlue);
ht->Draw();
hf->Draw("same");
li->Draw("same");
h2->Reset();
h1->Reset();

// bin 4
tfit->Draw("par[25]>>hf");
ttru->Draw("par[25]>>ht","step>300000");
ht->Scale(hf->Integral()/ht->Integral());
ht->SetLineColor(kBlue);
ht->Draw();
hf->Draw("same");
li->Draw("same");
h2->Reset();
h1->Reset();

// bin 5
tfit->Draw("par[25]>>hf");
ttru->Draw("par[25]>>ht","step>300000");
ht->Scale(hf->Integral()/ht->Integral());
ht->SetLineColor(kBlue);
ht->Draw();
hf->Draw("same");
li->Draw("same");
h2->Reset();
h1->Reset();

//h1->Scale(h2->Integral()/h1->Integral());
//h1->Draw();
//h2->Draw("same");
//li->Draw("same");

*/

}
