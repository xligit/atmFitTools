{

  TFile f("mcmctree.root");
  TFile g("likeplot.root");
 
  TGraph2D* hL = (TGraph2D*)g.Get("hL");
  TTree* path = (TTree*)f.Get("MCMCpath");

  hL->Draw("colz");
  path->Draw("par[1]:par[0]","","same");
}
