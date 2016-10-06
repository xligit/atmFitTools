{

TFile f("mcmc_steptest.root");
TFile g("demcmc_steptest.root");

TTree* mcmc = (TTree*)f.Get("MCMCpath");
TTree* demcmc = (TTree*)g.Get("MCMCpath");

TH2D* hlogm = new TH2D("hlogm","hlogm",1000,0,70000,5000,0,1600);
TH2D* hlogd = new TH2D("hlogd","hlogd",1000,0,70000,5000,0,1600);

mcmc->Draw("logL:step>>hlogm");
demcmc->Draw("logL:step>>hlogd");

hlogm->SetLineColor(kRed);
hlogd->SetLineColor(kBlue);
hlogm->SetMarkerColor(kRed);
hlogd->SetMarkerColor(kBlue);

hlogm->Draw();
hlogd->Draw("same");

}
