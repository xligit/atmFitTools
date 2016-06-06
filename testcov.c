{
// gROOT->ProcessLine(".L makeCov.cxx+");
 TFile f("mcmctree.root");
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov();
 maker->setParTree(tr);
 //gStyle->SetPalette(kBlackBody);
 maker->buildMatrix();
}
