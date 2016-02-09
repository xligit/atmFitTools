{
 gROOT->ProcessLine(".L makeCov.C+");
 TFile f("mcmctree.root");
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov();
 maker->setParTree(tr);
 maker->buildMatrix();
}
