{
 gROOT->ProcessLine(".L histoCompare.C++");
 gROOT->ProcessLine(".x ~/style.c");
 TChain atmmc("h1");
 atmmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/*stmusel_01*fQ.root");
 TChain atmdata("h1");
 atmdata.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run*.004*fQ.root");
 TH1F* hpdat = new TH1F("hpdat","hpdat",50,0,100);
 TH1F* hpmc = new TH1F("hpmc","hpmc",50,0,100);
 TH1F* hpiddat = new TH1F("hpiddat","hpiddat",50,-1000,600);
 TH1F* hpidmc = new TH1F("hpidmc","hpidmc",50,-1000,600);
 TH1F* hpiddatmu = new TH1F("hpiddatmu","hpiddatmu",50,-1000,600);
 TH1F* hpidmcmu = new TH1F("hpidmcmu","hpidmcmu",50,-1000,600);

// atmdata->Draw("fq1rmom[1][1]>>hpdat","(fqnse==2)&&(fqtotq[0]>100)&&(fq1rpcflg[1][1]==0)");
// atmmc->Draw("fq1rmom[1][1]>>hpmc","(fqnse==2)&&(fqtotq[0]>100)&&(fq1rpcflg[1][1]==0)");
 atmdata->Draw("fq1rnll[1][2]-fq1rnll[1][1]>>hpiddat","(fqnse==2)&&(fqtotq[0]>100)&&(fq1rpcflg[1][1]==0)");
 atmmc->Draw("fq1rnll[1][2]-fq1rnll[1][1]>>hpidmc","(fqnse==2)&&(fqtotq[0]>100)&&(fq1rpcflg[1][1]==0)");
 atmdata->Draw("fq1rnll[0][2]-fq1rnll[0][1]>>hpiddatmu","(fqnse==2)&&(fqtotq[0]>100)&&(fq1rpcflg[1][1]==0)");
 atmmc->Draw("fq1rnll[0][2]-fq1rnll[0][1]>>hpidmcmu","(fqnse==2)&&(fqtotq[0]>100)&&(fq1rpcflg[1][1]==0)");
// hpmc->Scale(hpdat->Integral()/hpmc->Integral());
 hpidmc->Scale(hpiddat->Integral()/hpidmc->Integral());
 hpidmcmu->Scale(hpiddatmu->Integral()/hpidmcmu->Integral());
 hpidmcmu->SetLineColor(kRed);
 hpidmc->SetLineColor(kBlue);

// hpidmc->Add(hpidmcmu);
 hpiddat->Add(hpiddatmu);

 histoCompare* hc = new histoCompare("comptest");
 hc->addHistogram(hpidmc,0);
 hc->addHistogram(hpidmcmu,0);
 hc->addHistogram(hpiddat,1);
// hc->addHistogram(hpidmcmu,0);
// hc->addHistogram(hpiddatmu,1);

}
