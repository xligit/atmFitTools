{
 gROOT->ProcessLine(".L makeCov.C+");
 TFile f("mcmc50000nom.root");
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov();
 maker->setParTree(tr);
 maker->buildMatrix();
}
