#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TGraph.h"


void wgtCompare()
{
  std::string j_path = "/disk2/sklb/t2kmc/14a/postfit_data_1p1h_biascorrection_20160310/weights/";
  std::string nu_path = "/disk2/usr5/xiaoyue/t2kmc14a/nu_mode/banffwgt/";
  std::string antinu_path = "/disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/banffwgt/";
  std::string nu_flavor[6] = {"numu_x_numu", "nue_x_nue", "numubar_x_numubar", "nuebar_x_nuebar", "numu_x_nue", "numubar_x_nuebar"};
  std::string antinu_flavor[6] = {"numu_crs_numu", "nue_crs_nue", "numubar_crs_numubar", "nuebar_crs_nuebar", "numu_crs_nue", "numubar_crs_nuebar"};
  float wgt1, wgt2;

  TFile *outfile = new TFile("wgtcompare.root", "recreate");

  // FHC
  for (int i = 0; i < 6; ++i) {
    TChain *c1 = new TChain("weightstree");
    TChain *c2 = new TChain("weightstree");
    c1->Add(Form("%s%s*weights.root", j_path.c_str(), nu_flavor[i].c_str()));
    c2->Add(Form("%s%s*weights.root", nu_path.c_str(), nu_flavor[i].c_str()));
    int nentries = c1->GetEntries();
    if (c2->GetEntries() != nentries) {
      std::cout<<nu_flavor[i]<<": old and new weightstrees don't have the same length!"<<std::endl;
      continue;
    }
    c1->SetBranchAddress("fWeight", &wgt1);
    c2->SetBranchAddress("fWeight", &wgt2);

    TH1D *hdiff = new TH1D(nu_flavor[i].c_str(), nu_flavor[i].c_str(), 200, -0.1, 0.1);
    TH2D *hdiffwgt1 = new TH2D(Form("%s_j",nu_flavor[i].c_str()), "", 200, 0, 2, 200, -0.1, 0.1);
    TH2D *hdiffwgt2 = new TH2D(Form("%s_x",nu_flavor[i].c_str()), "", 200, 0, 2, 200, -0.1, 0.1);
    TGraph *g12 = new TGraph(nentries);
    g12->SetName(Form("%s_g",nu_flavor[i].c_str()));

    for (int n = 0; n < nentries; ++n) {
      c1->GetEntry(n);
      c2->GetEntry(n);
      hdiff->Fill(wgt1-wgt2);
      hdiffwgt1->Fill(wgt1, wgt1-wgt2);
      hdiffwgt2->Fill(wgt2, wgt1-wgt2);
      g12->SetPoint(n, wgt1, wgt2);
    }

    outfile->cd();
    hdiff->Write();
    hdiffwgt1->Write();
    hdiffwgt2->Write();
    g12->Write();

    delete hdiff;
    delete hdiffwgt1;
    delete hdiffwgt2;
    delete g12;
  }

  // RHC
  for (int i = 0; i < 6; ++i) {
    TChain *c1 = new TChain("weightstree");
    TChain *c2 = new TChain("weightstree");
    c1->Add(Form("%s*%s*weights.root", j_path.c_str(), antinu_flavor[i].c_str()));
    c2->Add(Form("%s*%s*weights.root", antinu_path.c_str(), antinu_flavor[i].c_str()));
    int nentries = c1->GetEntries();
    if (c2->GetEntries() != nentries) {
      std::cout<<nu_flavor[i]<<": old and new weightstrees don't have the same length!"<<std::endl;
      continue;
    }
    c1->SetBranchAddress("fWeight", &wgt1);
    c2->SetBranchAddress("fWeight", &wgt2);

    TH1D *hdiff = new TH1D(antinu_flavor[i].c_str(), antinu_flavor[i].c_str(), 200, -0.1, 0.1);
    TH2D *hdiffwgt1 = new TH2D(Form("%s_j", antinu_flavor[i].c_str()), "", 200, 0, 2, 200, -0.1, 0.1);
    TH2D *hdiffwgt2 = new TH2D(Form("%s_x", antinu_flavor[i].c_str()), "", 200, 0, 2, 200, -0.1, 0.1);
    TGraph *g12 = new TGraph(nentries);
    g12->SetName(Form("%s_g", antinu_flavor[i].c_str()));
    for (int n = 0; n < nentries; ++n) {
      c1->GetEntry(n);
      c2->GetEntry(n);
      hdiff->Fill(wgt1-wgt2);
      hdiffwgt1->Fill(wgt1, wgt1-wgt2);
      hdiffwgt2->Fill(wgt2, wgt1-wgt2);
      g12->SetPoint(n, wgt1, wgt2);
    }

    outfile->cd();
    hdiff->Write();
    hdiffwgt1->Write();
    hdiffwgt2->Write();
    g12->Write();

    delete hdiff;
    delete hdiffwgt1;
    delete hdiffwgt2;
    delete g12;
  }

  outfile->Close();
}
