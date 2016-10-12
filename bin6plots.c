{
 TFile f("mcmctree.root");
 TTree* mcmc = (TTree*)f.Get("MCMCpath");

 TString rootname = "fit2";
 TString plotdir  = "~/transfer/";

 TString plotname = plotdir.Data();
 plotname.Append(rootname.Data());
 plotname.Append("_biaspar.png");
 mcmc->Draw("par[0]");
 c1->Print(plotname.Data());

/*
 plotname = plotdir.Data();
 plotname.Append(rootname.Data());
 plotname.Append("_biaspar_2d.png");
 mcmc->Draw("par[0]:par[1]","","colz");
 c1->Print(plotname.Data());

 plotname = plotdir.Data();
 plotname.Append(rootname.Data());
 plotname.Append("_biaspar_hpi_2d.png");
 mcmc->Draw("par[8]:par[1]","","colz");
 c1->Print(plotname.Data());
*/

}
