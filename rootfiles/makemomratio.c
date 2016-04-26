{
 TChain chm("h1");
 TChain chd("h1");
 chm.Add("cosmicmu_weight*ppmc*.root");
 chd.Add("cosmic*1s*ppdat*.root");
 TH1D* hpm = new TH1D("hpm","hpm",5000,0,10000);
 TH1D* hpd = new TH1D("hpd","hpd",5000,0,10000);
// chm.Draw("attribute[2]>>hpm","evtweight");
// chm.Draw("attribute[2]>>hpm");
// chd.Draw("attribute[2]>>hpd");
 hpm->Scale(hpd->Integral()/hpm->Integral());
 hpm->SetLineColor(kRed);
 hpm->SetLineWidth(3);
 hpd->SetLineWidth(3);
 hpm->SetStats(0);
 hpd->SetStats(0);
 hpm->Draw();
 hpd->Draw("same");

}
