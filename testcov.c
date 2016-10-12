{
 gROOT->ProcessLine(".L makeCov.cxx+");
// TFile f("./mcmc/closure_allpar_mcmcfit.root");
 TFile f("mcmctree.root");
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov("atmparsE.dat");
// makeCov *maker = new makeCov("fakepars1.dat");
 maker->setParTree(tr);
 maker->nburn = 40000;
 //gStyle->SetPalette(kBlackBody);
 maker->buildMatrix();
}
