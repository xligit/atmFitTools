#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TSpline.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObject.h"
#include "TAxis.h"
#include "TStopwatch.h"

#include "globalVars.h"

TH1D * nue_pid[2][9];
TH1D * nuebar_pid[2][9];
TH1D * numu_pid[2][9];
TH1D * numubar_pid[2][9];

void print_cov(TMatrixD*, TObjArray*);
int find_mode(int);
double xsec_wgt(int, TVectorD*);
double binned_xsec_wgt(int, int, int, int, bool, TVectorD*, int);
void draw_pid(TFile *);
int find_pid_bin(double, bool);
int find_etrue_bin(double ,bool);

#define np_maqe 21
#define np 13
#define nbins_pid_sub 60
#define nbins_pid_mult 90
#define nbins_etrue_sub 1
#define nbins_etrue_mult 1
TAxis *pid_sub = new TAxis(nbins_pid_sub, -3000, 3000);
TAxis *pid_mult = new TAxis(nbins_pid_mult, -3000, 6000);
// probably should change the binning of etrue;
TAxis *etrue_sub = new TAxis(nbins_etrue_sub, 0, 10);
TAxis *etrue_mult = new TAxis(nbins_etrue_mult, 1, 11);  
const int pidsubbins = nbins_pid_sub + 2;
const int pidmultbins = nbins_pid_mult + 2; 
const int etruesubbins = nbins_etrue_sub + 2;
const int etruemultbins = nbins_etrue_mult + 2;
std::map<int, int> modearray;

void MakeDistr(bool use_binned_spline = false, double sigma = 0, int nsub = 1) // nsub: number of sub-events, not to be confused with sub-GeV 
{
  TStopwatch *sw = new TStopwatch();
  gROOT->ProcessLine(".x ~/rootlogon.C");
  modearray[12]=0;
  modearray[-12]=1;
  modearray[14]=2;
  modearray[-14]=3;

  TH1::SetDefaultSumw2(true);
  
  TFile*infile = new TFile(filename.c_str(), "read");
  TTree *evtr = (TTree*)infile->Get("h1");
  TChain *sptr = new TChain("SplinesByEvent");
  sptr->Add(splfilename.c_str());
  
  evtr->SetBranchAddress("nring", &nring);
  evtr->SetBranchAddress("fqmrnring",fqmrnring);
  evtr->SetBranchAddress("fqnse", &fqnse);
  evtr->SetBranchAddress("nmue", &nmue);
  evtr->SetBranchAddress("nhitac", &nhitac);
  evtr->SetBranchAddress("ip", ip);
  evtr->SetBranchAddress("ipnu", ipnu);
  evtr->SetBranchAddress("mode", &mode);
  evtr->SetBranchAddress("evis", &evis);
  evtr->SetBranchAddress("dir", dir);
  evtr->SetBranchAddress("wall", &wall);
  evtr->SetBranchAddress("amomm", amomm);
  evtr->SetBranchAddress("amome", amome);
  evtr->SetBranchAddress("wallv", &wallv);
  evtr->SetBranchAddress("pnu", pnu);
  evtr->SetBranchAddress("fq1rmom", fq1rmom);
  evtr->SetBranchAddress("fq1rnll", fq1rnll);
  evtr->SetBranchAddress("fq1rdir", fq1rdir);
  evtr->SetBranchAddress("oscwgt", &oscwgt);
  evtr->SetBranchAddress("wgtosc", &wgtosc);
  evtr->SetBranchAddress("wgtflx", &wgtflx);
  evtr->SetBranchAddress("fqdwall", &fqdwall);
  evtr->SetBranchAddress("coszenith", &coszenith);
  evtr->SetBranchAddress("coslep", &coslep);
  sptr->SetBranchAddress("byEv_maqe_ccqe_gr", &byEv_maqe_ccqe_gr, &byEv_maqe_ccqe_br);
  sptr->SetBranchAddress("byEv_pfo_ccqe_gr", &byEv_pfo_ccqe_gr, &byEv_pfo_ccqe_br);
  sptr->SetBranchAddress("byEv_ebo_ccqe_gr", &byEv_ebo_ccqe_gr, &byEv_ebo_ccqe_br);
  sptr->SetBranchAddress("byEv_ca5_cc1pi_gr", &byEv_ca5_cc1pi_gr, &byEv_ca5_cc1pi_br);
  sptr->SetBranchAddress("byEv_ca5_ncpiz_gr", &byEv_ca5_ncpiz_gr, &byEv_ca5_ncpiz_br);
  sptr->SetBranchAddress("byEv_ca5_ncpipm_gr", &byEv_ca5_ncpipm_gr, &byEv_ca5_ncpipm_br);
  sptr->SetBranchAddress("byEv_manff_cc1pi_gr", &byEv_manff_cc1pi_gr, &byEv_manff_cc1pi_br);
  sptr->SetBranchAddress("byEv_manff_ncpiz_g", &byEv_manff_ncpiz_gr, &byEv_manff_ncpiz_br);
  sptr->SetBranchAddress("byEv_manff_ncpipm_gr", &byEv_manff_ncpipm_gr, &byEv_manff_ncpipm_br);
  sptr->SetBranchAddress("byEv_bgscl_cc1pi_gr", &byEv_bgscl_cc1pi_gr, &byEv_bgscl_cc1pi_br);
  sptr->SetBranchAddress("byEv_bgscl_ncpiz_gr", &byEv_bgscl_ncpiz_gr, &byEv_bgscl_ncpiz_br);
  sptr->SetBranchAddress("byEv_bgscl_ncpipm_gr", &byEv_bgscl_ncpipm_gr, &byEv_bgscl_ncpipm_br);
  sptr->SetBranchAddress("byEv_dismpishp_ccoth_gr", &byEv_dismpishp_ccoth_gr, &byEv_dismpishp_ccoth_br);
  sptr->SetBranchAddress("byEv_rpa_ccqe_gr", &byEv_rpa_ccqe_gr, &byEv_rpa_ccqe_br);     

  for (int i = 0; i < 9; ++i) {
    nue_pid[0][i] = new TH1D(Form("nue_%s_pid_sub",intnames[i].c_str()),"",60,-3000,3000);
    nue_pid[1][i] = new TH1D(Form("nue_%s_pid_mult",intnames[i].c_str()),"",90,-3000,6000);
    nuebar_pid[0][i] = new TH1D(Form("nuebar_%s_pid_sub",intnames[i].c_str()),"",60,-3000,3000);
    nuebar_pid[1][i] = new TH1D(Form("nuebar_%s_pid_mult",intnames[i].c_str()),"",90,-3000,6000);
    numu_pid[0][i] = new TH1D(Form("numu_%s_pid_sub",intnames[i].c_str()),"",60,-3000,3000);
    numu_pid[1][i] = new TH1D(Form("numu_%s_pid_mult",intnames[i].c_str()),"",90,-3000,6000);
    numubar_pid[0][i] = new TH1D(Form("numubar_%s_pid_sub",intnames[i].c_str()),"",60,-3000,3000);
    numubar_pid[1][i] = new TH1D(Form("numubar_%s_pid_mult",intnames[i].c_str()),"",90,-3000,6000);
  }  

  TFile *outfile = new TFile(Form("atm_distr_binned%d_%1.0fsigma_nsub%d.root",(int)use_binned_spline,sigma, nsub),"recreate");
  
  TFile *xsecf = new TFile(xsecfilename.c_str(), "read");
  TMatrixD *xsec_cov = (TMatrixD*)xsecf->Get("xsec_cov");  
  TObjArray *xsec_names = (TObjArray*)xsecf->Get("xsec_param_names");
  TVectorD *xsec_prior = (TVectorD*)xsecf->Get("xsec_param_prior");
  TVectorD *xsec_nom = (TVectorD*)xsecf->Get("xsec_param_nom");
  TVectorD *xsec_lb = (TVectorD*)xsecf->Get("xsec_param_lb");
  TVectorD *xsec_ub = (TVectorD*)xsecf->Get("xsec_param_ub");
  TVectorD *parvars = (TVectorD*)xsec_prior->Clone();  
  //print_cov(xsec_cov, xsec_names);
  double *sgm;
  sgm = new double[xsec_cov->GetNcols()];
  for (int i = 0; i < xsec_cov->GetNcols(); ++i) {
    sgm[i] = sqrt((*xsec_cov)(i,i));
    (*parvars)(i) += sigma * sgm[i] - (*xsec_nom)(i);
    if ((*parvars)(i) < (*xsec_lb)(i) - (*xsec_nom)(i)) (*parvars)(i) = (*xsec_lb)(i) - (*xsec_nom)(i);
    if ((*parvars)(i) > (*xsec_ub)(i) - (*xsec_nom)(i)) (*parvars)(i) = (*xsec_ub)(i) - (*xsec_nom)(i);
  }

  //xsec_prior->Print();
  //xsec_nom->Print();
  //parvars->Print();
  printf("%f(0), %f(4), %f(6), %f(7), %f(8), %f(9), %f(11)\n",(*parvars)(0),(*parvars)(4),(*parvars)(6),(*parvars)(7),(*parvars)(8),(*parvars)(9),(*parvars)(11));

  // --------------------------------------
  //   binned-splines reweight
  // --------------------------------------
  if (use_binned_spline) {
    nuespline = new TFile(Form("%s/nue_binned_splines_14a.root",splinepath.c_str()),"read");
    numuspline = new TFile(Form("%s/numu_binned_splines_14a.root",splinepath.c_str()),"read");
    nuebarspline = new TFile(Form("%s/nuebar_binned_splines_14a.root",splinepath.c_str()),"read");
    numubarspline = new TFile(Form("%s/numubar_binned_splines_14a.root",splinepath.c_str()),"read");

    Long64_t nentries = evtr->GetEntries();
    if (sptr->GetEntries() != nentries) {
      std::cerr<<"number of splines is not equal to number of event! quitting..."<<std::endl;
      exit(-1);
    }

    sw->Start();
    for (Long64_t n = 0; n < nentries; ++n) {
      evtr->GetEntry(n);
      //sptr->GetEntry(n);      
      if (nhitac>15 || evis<30 || wall < 200 || fqnse != nsub || fqmrnring[1]!=1) continue; // FCFV 1-ring 
      double wgt = wgtosc*wgtflx;
      int modee = find_mode(mode);
      switch (ipnu[0]) {
      case 12:
	if (evis<1000)    nue_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	else       nue_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	break;
      case -12:
	if (evis<1000)    nuebar_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	else      nuebar_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	break;
      case 14:
	if (evis<1000) 	  numu_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	else 	  numu_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	break;
      case -14:
	if (evis<1000) 	  numubar_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	else 	  numubar_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt);
	break;
      default:
	break;
      }      
    }
    sw->Print("m");

    sw->Stop(); sw->Reset(); sw->Start();
    for (int i = 0; i < nbins_pid_sub; ++i) {
      for (int j = 0; j < 9; ++j) {
	nue_pid[0][j]->SetBinContent(i+1, nue_pid[0][j]->GetBinContent(i+1)*binned_xsec_wgt(j,12,i+1,1,1,parvars,nsub));
	nue_pid[1][j]->SetBinContent(i+1, nue_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,12,i+1,1,0,parvars,nsub));
	nuebar_pid[0][j]->SetBinContent(i+1, nuebar_pid[0][j]->GetBinContent(i+1)*binned_xsec_wgt(j,-12,i+1,1,1,parvars,nsub));
	nuebar_pid[1][j]->SetBinContent(i+1, nuebar_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,-12,i+1,1,0,parvars,nsub));
	numu_pid[0][j]->SetBinContent(i+1, numu_pid[0][j]->GetBinContent(i+1)*binned_xsec_wgt(j,14,i+1,1,1,parvars,nsub));
	numu_pid[1][j]->SetBinContent(i+1, numu_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,14,i+1,1,0,parvars,nsub));
	numubar_pid[0][j]->SetBinContent(i+1, numubar_pid[0][j]->GetBinContent(i+1)*binned_xsec_wgt(j,-14,i+1,1,1,parvars,nsub));
	numubar_pid[1][j]->SetBinContent(i+1, numubar_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,-14,i+1,1,0,parvars,nsub));
	if (j==0) printf("%d, %f, %f, %f, %f\n", i, binned_xsec_wgt(j,12,i+1,1,1,parvars,nsub), binned_xsec_wgt(j,12,i+1,1,0,parvars,nsub), binned_xsec_wgt(j,-12,i+1,1,1,parvars,nsub), binned_xsec_wgt(j,-12,i+1,1,0,parvars,nsub));
      }
    }

    for (int i = nbins_pid_sub; i < nbins_pid_mult; ++i) {
      for (int j = 0; j < 9; ++j) {
	nue_pid[1][j]->SetBinContent(i+1, nue_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,12,i+1,1,0,parvars,nsub));
	nuebar_pid[1][j]->SetBinContent(i+1, nuebar_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,-12,i+1,1,0,parvars,nsub));
	numu_pid[1][j]->SetBinContent(i+1, numu_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,14,i+1,1,0,parvars,nsub));
	numubar_pid[1][j]->SetBinContent(i+1, numubar_pid[1][j]->GetBinContent(i+1)*binned_xsec_wgt(j,-14,i+1,1,0,parvars,nsub));
	if (j==0) printf("%d, %f, %f, %f, %f\n", i, binned_xsec_wgt(j,12,i+1,1,0,parvars,nsub), binned_xsec_wgt(j,12,i+1,1,0,parvars,nsub), binned_xsec_wgt(j,-12,i+1,1,0,parvars,nsub), binned_xsec_wgt(j,-12,i+1,1,0,parvars,nsub));
      }
    }
    sw->Print("m");

  }
  // ---------------------------------------  
  //    event-by-event reweight
  // ---------------------------------------
  else {
    std::string flv[] = {"nue","nuebar","numu","numubar"};
    std::string energy[] = {"subGeV", "multGeV"};
    for (int i = 0; i<4; ++i) {
	for (int k = 0; k<2; ++k) {
	  dmaqe[i][k] = new TH1D(Form("dmaqe_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s MaQE weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dpfo[i][k] = new TH1D(Form("dpfo_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s pF_O weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  debo[i][k] = new TH1D(Form("debo_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s Eb_O weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dca5_cc1pi[i][k] = new TH1D(Form("dca5_cc1pi_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s Ca5 cc1#pi weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dca5_ncpiz[i][k] = new TH1D(Form("dca5_ncpiz%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s Ca5 nc#pi^{0} weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dca5_ncpipm[i][k] = new TH1D(Form("dca5_ncpipm_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s Ca5 nc#pi^{#pm} weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dmanff_cc1pi[i][k] = new TH1D(Form("dmanff_cc1pi_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s MANFF cc1#pi weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dmanff_ncpiz[i][k] = new TH1D(Form("dmanff_ncpiz_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s MANFF nc#pi^{0} weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dmanff_ncpipm[i][k] = new TH1D(Form("dmanff_ncpipm_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s MANFF nc#pi^{#pm} weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dbgscl_cc1pi[i][k] = new TH1D(Form("dbgscl_cc1pi_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s BGSCL cc1#pi weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dbgscl_ncpiz[i][k] = new TH1D(Form("dbgscl_ncpiz_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s BGSCL nc#pi^{0} weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  dbgscl_ncpipm[i][k] = new TH1D(Form("dbgscl_ncpipm_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s BGSCL nc#pi^{#pm} weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	  ddismpishp[i][k] = new TH1D(Form("ddismpishp_%s_%dsubevt_%s",flv[i].c_str(),nsub,energy[k].c_str()),Form("%s DISMPIShp weight, %d sub-events, %s",flv[i].c_str(),nsub,energy[k].c_str()),100,0,2);
	}
    }

    Long64_t nentries = evtr->GetEntries();
    if (sptr->GetEntries() != nentries) {
      std::cerr<<"number of splines is not equal to number of event! quitting..."<<std::endl;
      exit(-1);
    }

    sw->Stop(); sw->Reset(); sw->Start();
    for (Long64_t n = 0; n < nentries; ++n) {
      evtr->GetEntry(n);
      sptr->GetEntry(n);
      
      if (nhitac>15 || evis<30 || wall < 200 || fqnse != nsub || fqmrnring[1]!=1) continue; // FCFV
      
      double wgt = wgtosc*wgtflx;
      int modee = find_mode(mode);
      switch (ipnu[0]) {
      case 12:
	if (evis<1000) nue_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	else nue_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	break;
      case -12:
	if (evis<1000) nuebar_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	else nuebar_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	break;
      case 14:
	if (evis<1000) numu_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	else numu_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	break;
      case -14:
	if (evis<1000) numubar_pid[0][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	else numubar_pid[1][modee]->Fill(fq1rnll[0][2]-fq1rnll[0][1], wgt*xsec_wgt(modee, parvars));
	break;
      }
      if (xsec_wgt(modee, parvars) < 0) printf("%d-th event has weight = %f < 0!\n", n, xsec_wgt(modee, parvars));
    }
  }
  sw->Print("m");
  
  for (int i = 0; i < 9; ++i) {
    nue_pid[0][i]->Scale(scale);
    nue_pid[1][i]->Scale(scale);
    nuebar_pid[0][i]->Scale(scale);
    nuebar_pid[1][i]->Scale(scale);
    numu_pid[0][i]->Scale(scale);
    numu_pid[1][i]->Scale(scale);
    numubar_pid[0][i]->Scale(scale);
    numubar_pid[1][i]->Scale(scale);
  }
  
  draw_pid(outfile);

  if (!use_binned_spline) {
    outfile->cd();
    for (int i = 0; i<4; ++i) {
      for (int k = 0; k<2; ++k) {
	dmaqe[i][k]->Write(0, TObject::kWriteDelete);
	dpfo[i][k]->Write(0, TObject::kWriteDelete);
	debo[i][k]->Write(0, TObject::kWriteDelete);
	dca5_cc1pi[i][k]->Write(0, TObject::kWriteDelete);
	dca5_ncpiz[i][k]->Write(0, TObject::kWriteDelete);
	dca5_ncpipm[i][k]->Write(0, TObject::kWriteDelete);
	dmanff_cc1pi[i][k]->Write(0, TObject::kWriteDelete);
	dmanff_ncpiz[i][k]->Write(0, TObject::kWriteDelete);
	dmanff_ncpipm[i][k]->Write(0, TObject::kWriteDelete);
	dbgscl_cc1pi[i][k]->Write(0, TObject::kWriteDelete);
	dbgscl_ncpiz[i][k]->Write(0, TObject::kWriteDelete);
	dbgscl_ncpipm[i][k]->Write(0, TObject::kWriteDelete);
	ddismpishp[i][k]->Write(0, TObject::kWriteDelete);
      }
    }
  }
  
  infile->Close();
  xsecf->Close();
  outfile->Close();

  if (nuespline!=NULL) delete nuespline;
  if (nuebarspline!=NULL) delete nuebarspline;
  if (numuspline!=NULL) delete numuspline;
  if (numubarspline!=NULL) delete numubarspline;

  /*
  for (int i = 0; i < 9; ++i) {
    if (nue_pid[0][i]!=NULL) delete nue_pid[0][i];
    if (nue_pid[1][i]!=NULL) delete nue_pid[1][i];
    if (nuebar_pid[0][i]!=NULL) delete nuebar_pid[0][i];
    if (nuebar_pid[1][i]!=NULL) delete nuebar_pid[1][i];
    if (numu_pid[0][i]!=NULL) delete numu_pid[0][i];
    if (numu_pid[1][i]!=NULL) delete numu_pid[1][i];
    if (numubar_pid[0][i]!=NULL) delete numubar_pid[0][i];
    if (numubar_pid[1][i]!=NULL) delete numubar_pid[1][i];
  }
  printf("deleted TH1Ds\n");
  */
  delete pid_sub;
  delete pid_mult;
  delete etrue_sub;
  delete etrue_mult;
  delete sw;
  if (sptr!=NULL) delete sptr;
  delete sgm;

}

int find_mode(int m)
{
  if (abs(m)==1) return 0;
  else if (abs(m)==2) return 8;
  else if (abs(m)>10 && abs(m)<14) return 1;
  else if (abs(m)==16) return 2;
  else if (abs(m)>16 && abs(m)<30) return 3;
  else if (abs(m)>30 && abs(m)<33) return 4;
  else if (abs(m)>32 && abs(m)<35) return 5;
  else if (abs(m)==36) return 6;
  else if (abs(m)>36 && abs(m)<53) return 7;
  else return -1;
}

double xsec_wgt(int m, TVectorD* par)
{
  double wgt = 1;
  switch (m) {
  case 0:
    wgt = byEv_maqe_ccqe_gr->Eval((*par)(0),0,"S")*byEv_pfo_ccqe_gr->Eval((*par)(4),0,"S")*byEv_ebo_ccqe_gr->Eval((*par)(6),0,"S");
    dmaqe[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_maqe_ccqe_gr->Eval((*par)(0),0,"S"));
    dpfo[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_pfo_ccqe_gr->Eval((*par)(4),0,"S"));
    debo[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_ebo_ccqe_gr->Eval((*par)(6),0,"S"));
    break;
  case 1:
    wgt = byEv_ca5_cc1pi_gr->Eval((*par)(7),0,"S")*byEv_manff_cc1pi_gr->Eval((*par)(8),0,"S")*byEv_bgscl_cc1pi_gr->Eval((*par)(9),0,"S");
    dca5_cc1pi[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_ca5_cc1pi_gr->Eval((*par)(7),0,"S"));
    dmanff_cc1pi[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_manff_cc1pi_gr->Eval((*par)(8),0,"S"));
    dbgscl_cc1pi[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_bgscl_cc1pi_gr->Eval((*par)(9),0,"S"));
    break;
  case 3:
    wgt = byEv_dismpishp_ccoth_gr->Eval((*par)(11),0,"S");
    ddismpishp[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(wgt);
    break;
  case 4:
    wgt = byEv_ca5_ncpiz_gr->Eval((*par)(7),0,"S")*byEv_manff_ncpiz_gr->Eval((*par)(8),0,"S")*byEv_bgscl_ncpiz_gr->Eval((*par)(9),0,"S");
    dca5_ncpiz[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_ca5_ncpiz_gr->Eval((*par)(7),0,"S"));
    dmanff_ncpiz[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_manff_ncpiz_gr->Eval((*par)(8),0,"S"));
    dbgscl_ncpiz[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_bgscl_ncpiz_gr->Eval((*par)(9),0,"S"));
    break;
  case 5:
    wgt = byEv_ca5_ncpipm_gr->Eval((*par)(7),0,"S")*byEv_manff_ncpipm_gr->Eval((*par)(8),0,"S")*byEv_bgscl_ncpipm_gr->Eval((*par)(9),0,"S");
    dca5_ncpipm[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_ca5_ncpipm_gr->Eval((*par)(7),0,"S"));
    dmanff_ncpipm[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_manff_ncpipm_gr->Eval((*par)(8),0,"S"));
    dbgscl_ncpipm[modearray[ipnu[0]]][(int)(evis<1000)]->Fill(byEv_bgscl_ncpipm_gr->Eval((*par)(9),0,"S"));
    break;
  default:
    break;
  }
  if (wgt<0) {
    //printf("xsec weight is < 0!\n");
    //exit(-1);
  }
  return wgt;
}

double binned_xsec_wgt(int m, int flavor, int i, int j, bool sub, TVectorD *par, int nsubev)
{
  double wgt = 1;
  TSpline *sp = NULL;
  TFile *file;
  if (flavor==12) file = nuespline;
  else if (flavor==-12) file = nuebarspline;
  else if (flavor==14) file = numuspline;
  else if (flavor==-14) file = numubarspline;
  else return 0;

  if (sub) { // sub-GeV
    switch (m) {
    case 0:
      sp = (TSpline3*)file->Get(Form("maqe_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("maqe_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(0));
      sp = (TSpline3*)file->Get(Form("pfo_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("pfo_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(4));
      sp = (TSpline3*)file->Get(Form("ebo_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ebo_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(6));
      break;
    case 1:
      sp = (TSpline3*)file->Get(Form("ca5_cc1pi_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ca5_cc1pi_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(7));
      sp = (TSpline3*)file->Get(Form("manff_cc1pi_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("manff_cc1pi_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(8));
      sp = (TSpline3*)file->Get(Form("bgscl_cc1pi_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("bgscl_cc1pi_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(9));
      break;
    case 3:
      sp = (TSpline3*)file->Get(Form("dismpishp_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("dismpishp_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(11));
      break;
    case 4:
      sp = (TSpline3*)file->Get(Form("ca5_ncpiz_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ca5_ncpiz_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(7));
      sp = (TSpline3*)file->Get(Form("manff_ncpiz_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("manff_ncpiz_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(8));
      sp = (TSpline3*)file->Get(Form("bgscl_ncpiz_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("bgscl_ncpiz_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(9));
      break;
    case 5:
      sp = (TSpline3*)file->Get(Form("ca5_ncpipm_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ca5_ncpipm_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(7));
      sp = (TSpline3*)file->Get(Form("manff_ncpipm_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("manff_ncpipm_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(8));
      sp = (TSpline3*)file->Get(Form("bgscl_ncpipm_p%d_e%d_ev%d_sub",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("bgscl_ncpipm_p%d_e%d_ev%d_sub",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(9));
      break;
    default:
      break;
    }
    
  } else { // multi-GeV
    switch (m) {
    case 0:
      sp = (TSpline3*)file->Get(Form("maqe_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("maqe_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(0));
      sp = (TSpline3*)file->Get(Form("pfo_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("pfo_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(4));
      sp = (TSpline3*)file->Get(Form("ebo_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ebo_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(6));
      break;
    case 1:
      sp = (TSpline3*)file->Get(Form("ca5_cc1pi_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ca5_cc1pi_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(7));
      sp = (TSpline3*)file->Get(Form("manff_cc1pi_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("manff_cc1pi_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(8));
      sp = (TSpline3*)file->Get(Form("bgscl_cc1pi_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("bgscl_cc1pi_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(9));
      break;
    case 3:
      sp = (TSpline3*)file->Get(Form("dismpishp_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("dismpishp_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(11));
      break;
    case 4:
      sp = (TSpline3*)file->Get(Form("ca5_ncpiz_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ca5_ncpiz_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(7));
      sp = (TSpline3*)file->Get(Form("manff_ncpiz_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("manff_ncpiz_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(8));
      sp = (TSpline3*)file->Get(Form("bgscl_ncpiz_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("bgscl_ncpiz_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(9));
      break;
    case 5:
      sp = (TSpline3*)file->Get(Form("ca5_ncpipm_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("ca5_ncpipm_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt = sp->Eval((*par)(7));
      sp = (TSpline3*)file->Get(Form("manff_ncpipm_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("manff_ncpipm_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(8));
      sp = (TSpline3*)file->Get(Form("bgscl_ncpipm_p%d_e%d_ev%d_mult",i,j,nsubev));
      if (sp==NULL) {std::cerr<<"cannot find TSpline3 "<<Form("bgscl_ncpipm_p%d_e%d_ev%d_mult",i,j,nsubev)<<std::endl; exit(-1); }
      wgt *= sp->Eval((*par)(9));
      break;
    default:
      break;
    }
  }
  if (wgt<0) {
    printf("binned xsec wgt = %f < 0!   %d flavor, %d-th bin, %d mode\n", wgt, flavor, i, m);
    //exit(-1);
  }
  return wgt;
}

void print_cov(TMatrixD *xsec_cov, TObjArray *xsec_names)
{
  static const int nxsecpars = xsec_cov->GetNcols();
  TString xsecnames[nxsecpars];
  for (int i = 0; i < nxsecpars; ++i) {
    xsecnames[i] = ((TObjString*)xsec_names->At(i))->GetString();
  }

  TH2D *h2 = new TH2D("h2","",nxsecpars,0,nxsecpars,nxsecpars,0,nxsecpars);
  for (int i = 0; i < nxsecpars; ++i) {
    for (int j = 0; j < nxsecpars; ++j) {
      h2->Fill(xsecnames[i].Data(), xsecnames[j].Data(), (*xsec_cov)(i, j));
    }
  }
  gStyle->SetOptStat(0);
  h2->GetXaxis()->SetNdivisions(200);
  h2->GetYaxis()->SetNdivisions(200);
  h2->GetXaxis()->LabelsOption("vv");
  h2->GetXaxis()->SetLabelSize(0.02);
  h2->GetYaxis()->SetLabelSize(0.02);  
  TCanvas *c = new TCanvas("c","",600,600);
  c->cd();
  h2->Draw("colz");

}

void draw_pid(TFile *f)
{
  TCanvas c1("c1","",600,600);
  TCanvas c2("c2","",600,600);
  // sub-GeV
  nue_pid[0][0]->SetName("nue_ccqe_sub");
  nue_pid[0][0]->SetFillColor(kBlue);
  nue_pid[0][0]->Add(nuebar_pid[0][0]);
  numu_pid[0][0]->SetName("numu_ccqe_sub");
  numu_pid[0][0]->SetFillColor(kRed);
  numu_pid[0][0]->Add(numubar_pid[0][0]);
  nue_pid[0][1]->SetName("nue_ccnqe_sub");
  nue_pid[0][1]->SetFillColor(kGreen);
  nue_pid[0][1]->Add(nue_pid[0][2]);
  nue_pid[0][1]->Add(nue_pid[0][3]);
  nue_pid[0][1]->Add(nue_pid[0][8]);
  nue_pid[0][1]->Add(nuebar_pid[0][1]);
  nue_pid[0][1]->Add(nuebar_pid[0][2]);
  nue_pid[0][1]->Add(nuebar_pid[0][3]);
  nue_pid[0][1]->Add(nuebar_pid[0][8]);
  numu_pid[0][1]->SetName("numu_ccnqe_sub");
  numu_pid[0][1]->SetFillColor(kYellow);
  numu_pid[0][1]->Add(numu_pid[0][2]);
  numu_pid[0][1]->Add(numu_pid[0][3]);
  numu_pid[0][1]->Add(numu_pid[0][8]);
  numu_pid[0][1]->Add(numubar_pid[0][1]);
  numu_pid[0][1]->Add(numubar_pid[0][2]);
  numu_pid[0][1]->Add(numubar_pid[0][3]);
  numu_pid[0][1]->Add(numubar_pid[0][8]);
  nue_pid[0][4]->SetName("nc_sub");
  nue_pid[0][4]->SetFillColor(kCyan);
  nue_pid[0][4]->Add(nue_pid[0][5]);
  nue_pid[0][4]->Add(nue_pid[0][6]);
  nue_pid[0][4]->Add(nue_pid[0][7]);
  nue_pid[0][4]->Add(numu_pid[0][4]);
  nue_pid[0][4]->Add(numu_pid[0][5]);
  nue_pid[0][4]->Add(numu_pid[0][6]);
  nue_pid[0][4]->Add(numu_pid[0][7]);
  nue_pid[0][4]->Add(nuebar_pid[0][4]);
  nue_pid[0][4]->Add(nuebar_pid[0][5]);
  nue_pid[0][4]->Add(nuebar_pid[0][6]);
  nue_pid[0][4]->Add(nuebar_pid[0][7]);
  nue_pid[0][4]->Add(numubar_pid[0][4]);
  nue_pid[0][4]->Add(numubar_pid[0][5]);
  nue_pid[0][4]->Add(numubar_pid[0][6]);
  nue_pid[0][4]->Add(numubar_pid[0][7]);
  nue_pid[0][4]->GetXaxis()->SetTitle("#Delta lnL");
  nue_pid[0][4]->SetNdivisions(505);
  THStack hs1("hs1","");
  hs1.SetHistogram(nue_pid[0][4]);
  hs1.Add(nue_pid[0][4]);
  hs1.Add(numu_pid[0][1]);
  hs1.Add(nue_pid[0][1]);
  hs1.Add(numu_pid[0][0]);
  hs1.Add(nue_pid[0][0]);
  TLegend lg1(0.15, 0.6, 0.4, 0.85);
  lg1.SetFillColor(0);
  lg1.AddEntry(nue_pid[0][0], "#nu_{e} CCQE", "f");
  lg1.AddEntry(numu_pid[0][0], "#nu_{#mu} CCQE", "f");
  lg1.AddEntry(nue_pid[0][1], "#nu_{e} CCnQE", "f");
  lg1.AddEntry(numu_pid[0][1], "#nu_{#mu} CCnQE", "f");
  lg1.AddEntry(nue_pid[0][4], "NC", "f");
  c1.cd();
  hs1.Draw();
  lg1.Draw("same");
  // multi-GeV
  nue_pid[1][0]->SetName("nue_ccqe_mult");
  nue_pid[1][0]->SetFillColor(kBlue);
  nue_pid[1][0]->Add(nuebar_pid[1][0]);
  numu_pid[1][0]->SetName("numu_ccqe_mult");
  numu_pid[1][0]->SetFillColor(kRed);
  numu_pid[1][0]->Add(numubar_pid[1][0]);
  nue_pid[1][1]->SetName("nue_ccnqe_mult");
  nue_pid[1][1]->SetFillColor(kGreen);
  nue_pid[1][1]->Add(nue_pid[1][2]);
  nue_pid[1][1]->Add(nue_pid[1][3]);
  nue_pid[1][1]->Add(nue_pid[1][8]);
  nue_pid[1][1]->Add(nuebar_pid[1][1]);
  nue_pid[1][1]->Add(nuebar_pid[1][2]);
  nue_pid[1][1]->Add(nuebar_pid[1][3]);
  nue_pid[1][1]->Add(nuebar_pid[1][8]);
  numu_pid[1][1]->SetName("numu_ccnqe_mult");
  numu_pid[1][1]->SetFillColor(kYellow);
  numu_pid[1][1]->Add(numu_pid[1][2]);
  numu_pid[1][1]->Add(numu_pid[1][3]);
  numu_pid[1][1]->Add(numu_pid[1][8]);
  numu_pid[1][1]->Add(numubar_pid[1][1]);
  numu_pid[1][1]->Add(numubar_pid[1][2]);
  numu_pid[1][1]->Add(numubar_pid[1][3]);
  numu_pid[1][1]->Add(numubar_pid[1][8]);
  nue_pid[1][4]->SetName("nc_mult");
  nue_pid[1][4]->SetFillColor(kCyan);
  nue_pid[1][4]->Add(nue_pid[1][5]);
  nue_pid[1][4]->Add(nue_pid[1][6]);
  nue_pid[1][4]->Add(nue_pid[1][7]);
  nue_pid[1][4]->Add(numu_pid[1][4]);
  nue_pid[1][4]->Add(numu_pid[1][5]);
  nue_pid[1][4]->Add(numu_pid[1][6]);
  nue_pid[1][4]->Add(numu_pid[1][7]);
  nue_pid[1][4]->Add(nuebar_pid[1][4]);
  nue_pid[1][4]->Add(nuebar_pid[1][5]);
  nue_pid[1][4]->Add(nuebar_pid[1][6]);
  nue_pid[1][4]->Add(nuebar_pid[1][7]);
  nue_pid[1][4]->Add(numubar_pid[1][4]);
  nue_pid[1][4]->Add(numubar_pid[1][5]);
  nue_pid[1][4]->Add(numubar_pid[1][6]);
  nue_pid[1][4]->Add(numubar_pid[1][7]);
  nue_pid[1][4]->GetXaxis()->SetTitle("#Delta lnL");
  nue_pid[1][4]->SetNdivisions(505);
  THStack hs2("hs2","");
  hs2.SetHistogram(nue_pid[1][4]);
  hs2.Add(nue_pid[1][4]);
  hs2.Add(numu_pid[1][1]);
  hs2.Add(nue_pid[1][1]);
  hs2.Add(numu_pid[1][0]);
  hs2.Add(nue_pid[1][0]);
  TLegend lg2(0.6, 0.6, 0.85, 0.85);
  lg2.SetFillColor(0);
  lg2.AddEntry(nue_pid[1][0], "#nu_{e} CCQE", "f");
  lg2.AddEntry(numu_pid[1][0], "#nu_{#mu} CCQE", "f");
  lg2.AddEntry(nue_pid[1][1], "#nu_{e} CCnQE", "f");
  lg2.AddEntry(numu_pid[1][1], "#nu_{#mu} CCnQE", "f");
  lg2.AddEntry(nue_pid[1][4], "NC", "f");
  c2.cd();
  hs2.Draw();
  lg2.Draw("same");
 
  f->cd();
  c1.Write(0, TObject::kWriteDelete);
  c2.Write(0, TObject::kWriteDelete);
  nue_pid[0][0]->Write(0, TObject::kWriteDelete);
  numu_pid[0][0]->Write(0, TObject::kWriteDelete);
  nue_pid[0][1]->Write(0, TObject::kWriteDelete);
  numu_pid[0][1]->Write(0, TObject::kWriteDelete);
  nue_pid[0][4]->Write(0, TObject::kWriteDelete);
  nue_pid[1][0]->Write(0, TObject::kWriteDelete);
  numu_pid[1][0]->Write(0, TObject::kWriteDelete);
  nue_pid[1][1]->Write(0, TObject::kWriteDelete);
  numu_pid[1][1]->Write(0, TObject::kWriteDelete);
  nue_pid[1][4]->Write(0, TObject::kWriteDelete);

}

int find_pid_bin(double pidpar, bool sub)
{
  if (sub) return pid_sub->FindBin(pidpar);
  else return pid_mult->FindBin(pidpar);
}

int find_etrue_bin(double etrue, bool sub)
{
  if (sub) return etrue_sub->FindBin(etrue);
  else return etrue_mult->FindBin(etrue);
}
