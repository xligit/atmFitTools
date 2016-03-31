#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
#include "TAxis.h"
#include "TH2D.h"

#include "globalVars.h"

#define np_maqe 21
#define np 13
#define nbins_pid_sub 60
#define nbins_pid_mult 90
#define nbins_etrue_sub 1
#define nbins_etrue_mult 1

int find_mode(int);
int find_pid_bin(double, bool);
int find_etrue_bin(double, bool);

TAxis *pid_sub = new TAxis(nbins_pid_sub, -3000, 3000);
TAxis *pid_mult = new TAxis(nbins_pid_mult, -3000, 6000);
// probably should change the binning of etrue;
TAxis *etrue_sub = new TAxis(nbins_etrue_sub, 0, 10);
TAxis *etrue_mult = new TAxis(nbins_etrue_mult, 1, 11);  
const int pidsubbins = nbins_pid_sub + 2;
const int pidmultbins = nbins_pid_mult + 2; 
const int etruesubbins = nbins_etrue_sub + 2;
const int etruemultbins = nbins_etrue_mult + 2;

void MakeBinnedSplines(const std::string flavor = "", int subev = 0)
{
  
  gROOT->ProcessLine(".x ~/rootlogon.C");

  // Initialize
  TFile *outfile;
  TFile *plotfile;
  int nucode;
  if (flavor=="nue") {
    outfile = new TFile(Form("%s/nue_binned_splines_14a.root",splinepath.c_str()),"recreate");
    plotfile = new TFile(Form("%s/nue_binned_splines_%dsubev_plot_14a.root",splinepath.c_str(),subev+1),"recreate");
    nucode = 12;
  } else if (flavor=="numu") {
    outfile = new TFile(Form("%s/numu_binned_splines_14a.root",splinepath.c_str()),"recreate");
    plotfile = new TFile(Form("%s/numu_binned_splines_%dsubev_plot_14a.root",splinepath.c_str(),subev+1),"recreate");   
    nucode = 14;
  } else if (flavor=="nuebar") {
    outfile = new TFile(Form("%s/nuebar_binned_splines_14a.root",splinepath.c_str()),"recreate");
    plotfile = new TFile(Form("%s/nuebar_binned_splines_%dsubev_plot_14a.root",splinepath.c_str(),subev+1),"recreate");
    nucode = -12;
  } else if (flavor=="numubar") {
    outfile = new TFile(Form("%s/numubar_binned_splines_14a.root",splinepath.c_str()),"recreate");
    plotfile = new TFile(Form("%s/numubar_binned_splines_%dsubev_plot_14a.root",splinepath.c_str(),subev+1),"recreate");
    nucode = -14;
  } else {
    printf("Wrong flavor input. Please enter one from the following:\n");
    printf("nue, nuebar, numu, numbar.\n");
    exit(-1);
  }

  TFile*infile = new TFile(filename.c_str(), "read");
  TTree *evtr = (TTree*)infile->Get("h1");
  TChain *sptr = new TChain("SplinesByEvent");
  sptr->Add(splfilename.c_str());
  
  evtr->SetBranchAddress("nring", &nring);
  evtr->SetBranchAddress("fqmrnring",fqmrnring);
  evtr->SetBranchAddress("fqnse", &fqnse);
  evtr->SetBranchAddress("nmue", &nmue);
  evtr->SetBranchAddress("nhitac", &nhitac);
  //evtr->SetBranchAddress("evclass", &evclass);
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
  
  
  Long64_t nentries = evtr->GetEntries();
  if (sptr->GetEntries() != nentries) {
    std::cerr<<"number of splines is not equal to number of event! quitting..."<<std::endl;
    exit(-1);
  }

  // prepare output

  TSpline3 *g_maqe_sub[pidsubbins][etruesubbins][2];// one sub-event and two sub-events
  TSpline3 *g_maqe_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_pfo_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_pfo_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_ebo_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_ebo_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_ca5_cc1pi_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_ca5_cc1pi_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_ca5_ncpiz_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_ca5_ncpiz_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_ca5_ncpipm_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_ca5_ncpipm_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_manff_cc1pi_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_manff_cc1pi_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_manff_ncpiz_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_manff_ncpiz_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_manff_ncpipm_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_manff_ncpipm_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_bgscl_cc1pi_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_bgscl_cc1pi_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_bgscl_ncpiz_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_bgscl_ncpiz_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_bgscl_ncpipm_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_bgscl_ncpipm_mult[pidmultbins][etruemultbins][2];
  TSpline3 *g_dismpishp_sub[pidsubbins][etruesubbins][2];
  TSpline3 *g_dismpishp_mult[pidmultbins][etruemultbins][2];

  double maqe_sub[pidsubbins][etruesubbins][2][np_maqe] = {0.};// one sub-event and two sub-events
  double maqe_mult[pidmultbins][etruemultbins][2][np_maqe] = {0.};
  double pfo_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double pfo_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double ebo_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double ebo_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double ca5_cc1pi_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double ca5_cc1pi_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double ca5_ncpiz_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double ca5_ncpiz_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double ca5_ncpipm_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double ca5_ncpipm_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double manff_cc1pi_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double manff_cc1pi_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double manff_ncpiz_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double manff_ncpiz_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double manff_ncpipm_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double manff_ncpipm_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double bgscl_cc1pi_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double bgscl_cc1pi_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double bgscl_ncpiz_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double bgscl_ncpiz_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double bgscl_ncpipm_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double bgscl_ncpipm_mult[pidmultbins][etruemultbins][2][np] = {0.};
  double dismpishp_sub[pidsubbins][etruesubbins][2][np] = {0.};
  double dismpishp_mult[pidmultbins][etruemultbins][2][np] = {0.};

  double ccqe_sub[pidsubbins][etruesubbins][2] = {0.};// one sub-event and two sub-events
  double ccqe_mult[pidmultbins][etruemultbins][2] = {0.};
  double cc1pi_sub[pidsubbins][etruesubbins][2] = {0.};
  double cc1pi_mult[pidmultbins][etruemultbins][2] = {0.};
  double ncpiz_sub[pidsubbins][etruesubbins][2] = {0.};
  double ncpiz_mult[pidmultbins][etruemultbins][2] = {0.};
  double ncpipm_sub[pidsubbins][etruesubbins][2] = {0.};
  double ncpipm_mult[pidmultbins][etruemultbins][2] = {0.};
  double ccoth_sub[pidsubbins][etruesubbins][2] = {0.};
  double ccoth_mult[pidmultbins][etruemultbins][2] = {0.};
  
  double maqe[np_maqe], pfo[np], ebo[np], ca5_cc1pi[np], ca5_ncpiz[np], ca5_ncpipm[np], manff_cc1pi[np], manff_ncpiz[np], manff_ncpipm[np], bgscl_cc1pi[np], bgscl_ncpiz[np], bgscl_ncpipm[np], dismpishp[np];
  sptr->GetEntry(0);
  double xy[2];
  for (int i = 0; i < np; ++i) {
    byEv_maqe_ccqe_gr->GetPoint(i,xy[0],xy[1]);
    maqe[i] = xy[0];
    byEv_pfo_ccqe_gr->GetPoint(i,xy[0],xy[1]);
    pfo[i] = xy[0];
    byEv_ebo_ccqe_gr->GetPoint(i,xy[0],xy[1]);
    ebo[i] = xy[0];
    byEv_ca5_cc1pi_gr->GetPoint(i,xy[0],xy[1]);
    ca5_cc1pi[i] = xy[0];
    byEv_ca5_ncpiz_gr->GetPoint(i,xy[0],xy[1]);
    ca5_ncpiz[i] = xy[0];
    byEv_ca5_ncpipm_gr->GetPoint(i,xy[0],xy[1]);
    ca5_ncpipm[i] = xy[0];
    byEv_manff_cc1pi_gr->GetPoint(i,xy[0],xy[1]);
    manff_cc1pi[i] = xy[0];
    byEv_manff_ncpiz_gr->GetPoint(i,xy[0],xy[1]);
    manff_ncpiz[i] = xy[0];
    byEv_manff_ncpipm_gr->GetPoint(i,xy[0],xy[1]);
    manff_ncpipm[i] = xy[0];
    byEv_bgscl_cc1pi_gr->GetPoint(i,xy[0],xy[1]);
    bgscl_cc1pi[i] = xy[0];
    byEv_bgscl_ncpiz_gr->GetPoint(i,xy[0],xy[1]);
    bgscl_ncpiz[i] = xy[0];
    byEv_bgscl_ncpipm_gr->GetPoint(i,xy[0],xy[1]);
    bgscl_ncpipm[i] = xy[0];
    byEv_dismpishp_ccoth_gr->GetPoint(i,xy[0],xy[1]);
    dismpishp[i] = xy[0];    
  }
  for (int i = np; i < np_maqe; ++i) {
    byEv_maqe_ccqe_gr->GetPoint(i,xy[0],xy[1]);
    maqe[i] = xy[0];
  }
  
  for (Long64_t n = 0; n < nentries; ++n) {
    evtr->GetEntry(n);
    sptr->GetEntry(n);
    if (nhitac>15 || evis<30 || wall<200 || ipnu[0]!=nucode || fqnse>2 || fqmrnring[1]!=1) continue;
    double wgt = wgtosc*wgtflx;
    int modee = find_mode(mode);
    int ebin = fqnse - 1;
    int ibin, jbin;
    double x, y;
    
    if (evis<1000) { // sub-GeV
      ibin = find_pid_bin(fq1rnll[0][2]-fq1rnll[0][1], 1);
      jbin = 1;//find_etrue_bin(pnu[0], 1);
      switch (modee) {
      case 0:
	for (int k = 0; k < np; ++k) {
	  byEv_maqe_ccqe_gr->GetPoint(k,x,y);
	  maqe_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_pfo_ccqe_gr->GetPoint(k,x,y);
	  pfo_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_ebo_ccqe_gr->GetPoint(k,x,y);
	  ebo_sub[ibin][jbin][ebin][k] += y*wgt;
	}
	for (int k = np; k < np_maqe; ++k) {
	  byEv_maqe_ccqe_gr->GetPoint(k,x,y);
	  maqe_sub[ibin][jbin][ebin][k] += y*wgt;
	}
	ccqe_sub[ibin][jbin][ebin] += wgt;
	break;
      case 1:
	for (int k = 0; k < np; ++k) {
	  byEv_ca5_cc1pi_gr->GetPoint(k,x,y);
	  ca5_cc1pi_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_manff_cc1pi_gr->GetPoint(k,x,y);
	  manff_cc1pi_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_bgscl_cc1pi_gr->GetPoint(k,x,y);
	  bgscl_cc1pi_sub[ibin][jbin][ebin][k] += y*wgt;
	}
	cc1pi_sub[ibin][jbin][ebin] += wgt;
	break;
      case 4:
	for (int k = 0; k < np; ++k) {
	  byEv_ca5_ncpiz_gr->GetPoint(k,x,y);
	  ca5_ncpiz_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_manff_ncpiz_gr->GetPoint(k,x,y);
	  manff_ncpiz_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_bgscl_ncpiz_gr->GetPoint(k,x,y);
	  bgscl_ncpiz_sub[ibin][jbin][ebin][k] += y*wgt;
	}
	ncpiz_sub[ibin][jbin][ebin] += wgt;
	break;
      case 5:
	for (int k = 0; k < np; ++k) {
	  byEv_ca5_ncpipm_gr->GetPoint(k,x,y);
	  ca5_ncpipm_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_manff_ncpipm_gr->GetPoint(k,x,y);
	  manff_ncpipm_sub[ibin][jbin][ebin][k] += y*wgt;
	  byEv_bgscl_ncpipm_gr->GetPoint(k,x,y);
	  bgscl_ncpipm_sub[ibin][jbin][ebin][k] += y*wgt;
	}
	ncpipm_sub[ibin][jbin][ebin] += wgt;
	break;
      case 3:
	for (int k = 0; k < np; ++k) {
	  byEv_dismpishp_ccoth_gr->GetPoint(k,x,y);
	  dismpishp_sub[ibin][jbin][ebin][k] += y*wgt;
	}
	ccoth_sub[ibin][jbin][ebin] += wgt;
	break;
      default:
	break;
      }
      
    } else { // multi-GeV
      ibin = find_pid_bin(fq1rnll[0][2]-fq1rnll[0][1], 0);
      jbin = 1;//find_etrue_bin(pnu[0], 0);
      switch (modee) {
      case 0:
	for (int k = 0; k < np; ++k) {
	  byEv_maqe_ccqe_gr->GetPoint(k,x,y);
	  maqe_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_pfo_ccqe_gr->GetPoint(k,x,y);
	  pfo_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_ebo_ccqe_gr->GetPoint(k,x,y);
	  ebo_mult[ibin][jbin][ebin][k] += y*wgt;
	}
	for (int k = np; k < np_maqe; ++k) {
	  byEv_maqe_ccqe_gr->GetPoint(k,x,y);
	  maqe_mult[ibin][jbin][ebin][k] += y*wgt;
	}
	ccqe_mult[ibin][jbin][ebin] += wgt;
	break;
      case 1:
	for (int k = 0; k < np; ++k) {
	  byEv_ca5_cc1pi_gr->GetPoint(k,x,y);
	  ca5_cc1pi_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_manff_cc1pi_gr->GetPoint(k,x,y);
	  manff_cc1pi_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_bgscl_cc1pi_gr->GetPoint(k,x,y);
	  bgscl_cc1pi_mult[ibin][jbin][ebin][k] += y*wgt;
	}
	cc1pi_mult[ibin][jbin][ebin] += wgt;
	break;
      case 4:
	for (int k = 0; k < np; ++k) {
	  byEv_ca5_ncpiz_gr->GetPoint(k,x,y);
	  ca5_ncpiz_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_manff_ncpiz_gr->GetPoint(k,x,y);
	  manff_ncpiz_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_bgscl_ncpiz_gr->GetPoint(k,x,y);
	  bgscl_ncpiz_mult[ibin][jbin][ebin][k] += y*wgt;
	}
	ncpiz_mult[ibin][jbin][ebin] += wgt;
	break;
      case 5:
	for (int k = 0; k < np; ++k) {
	  byEv_ca5_ncpipm_gr->GetPoint(k,x,y);
	  ca5_ncpipm_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_manff_ncpipm_gr->GetPoint(k,x,y);
	  manff_ncpipm_mult[ibin][jbin][ebin][k] += y*wgt;
	  byEv_bgscl_ncpipm_gr->GetPoint(k,x,y);
	  bgscl_ncpipm_mult[ibin][jbin][ebin][k] += y*wgt;
	}
	ncpipm_mult[ibin][jbin][ebin] += wgt;
	break;
      case 3:
	for (int k = 0; k < np; ++k) {
	  byEv_dismpishp_ccoth_gr->GetPoint(k,x,y);
	  dismpishp_mult[ibin][jbin][ebin][k] += y*wgt;
	}
	ccoth_mult[ibin][jbin][ebin] += wgt;
	break;
      default:
	break;
      }      
    }    
  } // end of event loop

  // TH2Ds to illustrate binned splines
  TH2D *h_maqe_sub = new TH2D("h_maqe_sub","MAQE sub-GeV", np_maqe, (3*maqe[0]-maqe[1])/2., (3*maqe[20]-maqe[19])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_maqe_mult = new TH2D("h_maqe_mult","MAQE mult-GeV", np_maqe, (3*maqe[0]-maqe[1])/2., (3*maqe[20]-maqe[19])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_pfo_sub = new TH2D("h_pfo_sub","pfO sub-GeV", np, (3*pfo[0]-pfo[1])/2., (3*pfo[12]-pfo[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_pfo_mult = new TH2D("h_pfo_mult","pfO mult-GeV", np, (3*pfo[0]-pfo[1])/2., (3*pfo[12]-pfo[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_ebo_sub = new TH2D("h_ebo_sub","EbO sub-GeV", np, (3*ebo[0]-ebo[1])/2., (3*ebo[12]-ebo[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_ebo_mult = new TH2D("h_ebo_mult","EbO mult-GeV", np, (3*ebo[0]-ebo[1])/2., (3*ebo[12]-ebo[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_ca5_cc1pi_sub = new TH2D("h_ca5_cc1pi_sub","CA5_CC1PI sub-GeV", np, (3*ca5_cc1pi[0]-ca5_cc1pi[1])/2., (3*ca5_cc1pi[12]-ca5_cc1pi[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_ca5_cc1pi_mult = new TH2D("h_ca5_cc1pi_mult","CA5_CC1PI mult-GeV", np, (3*ca5_cc1pi[0]-ca5_cc1pi[1])/2., (3*ca5_cc1pi[12]-ca5_cc1pi[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_ca5_ncpiz_sub = new TH2D("h_ca5_ncpiz_sub","CA5_NCPIZ sub-GeV", np, (3*ca5_ncpiz[0]-ca5_ncpiz[1])/2., (3*ca5_ncpiz[12]-ca5_ncpiz[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_ca5_ncpiz_mult = new TH2D("h_ca5_ncpiz_mult","CA5_NCPIZ mult-GeV", np, (3*ca5_ncpiz[0]-ca5_ncpiz[1])/2., (3*ca5_ncpiz[12]-ca5_ncpiz[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_ca5_ncpipm_sub = new TH2D("h_ca5_ncpipm_sub","CA5_NCPIPM sub-GeV", np, (3*ca5_ncpipm[0]-ca5_ncpipm[1])/2., (3*ca5_ncpipm[12]-ca5_ncpipm[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_ca5_ncpipm_mult = new TH2D("h_ca5_ncpipm_mult","CA5_NCPIPM mult-GeV", np, (3*ca5_ncpipm[0]-ca5_ncpipm[1])/2., (3*ca5_ncpipm[12]-ca5_ncpipm[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_manff_cc1pi_sub = new TH2D("h_manff_cc1pi_sub","MANFF_CC1PI sub-GeV", np, (3*manff_cc1pi[0]-manff_cc1pi[1])/2., (3*manff_cc1pi[12]-manff_cc1pi[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_manff_cc1pi_mult = new TH2D("h_manff_cc1pi_mult","MANFF_CC1PI mult-GeV", np, (3*manff_cc1pi[0]-manff_cc1pi[1])/2., (3*manff_cc1pi[12]-manff_cc1pi[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_manff_ncpiz_sub = new TH2D("h_manff_ncpiz_sub","MANFF_NCPIZ sub-GeV", np, (3*manff_ncpiz[0]-manff_ncpiz[1])/2., (3*manff_ncpiz[12]-manff_ncpiz[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_manff_ncpiz_mult = new TH2D("h_manff_ncpiz_mult","MANFF_NCPIZ mult-GeV", np, (3*manff_ncpiz[0]-manff_ncpiz[1])/2., (3*manff_ncpiz[12]-manff_ncpiz[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_manff_ncpipm_sub = new TH2D("h_manff_ncpipm_sub","MANFF_NCPIPM sub-GeV", np, (3*manff_ncpipm[0]-manff_ncpipm[1])/2., (3*manff_ncpipm[12]-manff_ncpipm[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_manff_ncpipm_mult = new TH2D("h_manff_ncpipm_mult","MANFF_NCPIPM mult-GeV", np, (3*manff_ncpipm[0]-manff_ncpipm[1])/2., (3*manff_ncpipm[12]-manff_ncpipm[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_bgscl_cc1pi_sub = new TH2D("h_bgscl_cc1pi_sub","BGSCL_CC1PI sub-GeV", np, (3*bgscl_cc1pi[0]-bgscl_cc1pi[1])/2., (3*bgscl_cc1pi[12]-bgscl_cc1pi[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_bgscl_cc1pi_mult = new TH2D("h_bgscl_cc1pi_mult","BGSCL_CC1PI mult-GeV", np, (3*bgscl_cc1pi[0]-bgscl_cc1pi[1])/2., (3*bgscl_cc1pi[12]-bgscl_cc1pi[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_bgscl_ncpiz_sub = new TH2D("h_bgscl_ncpiz_sub","BGSCL_NCPIZ sub-GeV", np, (3*bgscl_ncpiz[0]-bgscl_ncpiz[1])/2., (3*bgscl_ncpiz[12]-bgscl_ncpiz[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_bgscl_ncpiz_mult = new TH2D("h_bgscl_ncpiz_mult","BGSCL_NCPIZ mult-GeV", np, (3*bgscl_ncpiz[0]-bgscl_ncpiz[1])/2., (3*bgscl_ncpiz[12]-bgscl_ncpiz[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_bgscl_ncpipm_sub = new TH2D("h_bgscl_ncpipm_sub","BGSCL_NCPIPM sub-GeV", np, (3*bgscl_ncpipm[0]-bgscl_ncpipm[1])/2., (3*bgscl_ncpipm[12]-bgscl_ncpipm[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_bgscl_ncpipm_mult = new TH2D("h_bgscl_ncpipm_mult","BGSCL_NCPIPM mult-GeV", np, (3*bgscl_ncpipm[0]-bgscl_ncpipm[1])/2., (3*bgscl_ncpipm[12]-bgscl_ncpipm[11])/2., nbins_pid_mult, -3000, 6000);
  TH2D *h_dismpishp_sub = new TH2D("h_dismpishp_sub","DISMPISHP sub-GeV", np, (3*dismpishp[0]-dismpishp[1])/2., (3*dismpishp[12]-dismpishp[11])/2., nbins_pid_sub, -3000, 3000);
  TH2D *h_dismpishp_mult = new TH2D("h_dismpishp_mult","DISMPISHP mult-GeV", np, (3*dismpishp[0]-dismpishp[1])/2., (3*dismpishp[12]-dismpishp[11])/2., nbins_pid_mult, -3000, 6000);
  
  // construct and save splines to file -- sub-GeV  
  for (int i = 0; i < pidsubbins; ++i) {
    for (int j = 1; j < 2/*etruesubbins*/; ++j) {
      for (int k = 0; k < 2; ++k) {

	if (ccqe_sub[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    maqe_sub[i][j][k][e] /= ccqe_sub[i][j][k];
	    pfo_sub[i][j][k][e] /= ccqe_sub[i][j][k];
	    ebo_sub[i][j][k][e] /= ccqe_sub[i][j][k];
	  }
	  for (int e = np; e < np_maqe; ++e) {
	    maqe_sub[i][j][k][e] /= ccqe_sub[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    maqe_sub[i][j][k][e] = 1;
	    pfo_sub[i][j][k][e] = 1;
	    ebo_sub[i][j][k][e] = 1;
	  }
	  for (int e = np; e < np_maqe; ++e) {
	    maqe_sub[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  maqe_sub[i][j][k][e] /= maqe_sub[i][j][k][14];
	  pfo_sub[i][j][k][e] /= pfo_sub[i][j][k][6];
	  ebo_sub[i][j][k][e] /= ebo_sub[i][j][k][6];
	  if (k==subev) {
	    h_maqe_sub->SetBinContent(e+1, i+1, maqe_sub[i][j][k][e]);
	    h_pfo_sub->SetBinContent(e+1, i+1, pfo_sub[i][j][k][e]);
	    h_ebo_sub->SetBinContent(e+1, i+1, ebo_sub[i][j][k][e]);
	  }
	}
	for (int e = np; e < np_maqe; ++e) {
	  maqe_sub[i][j][k][e] /=maqe_sub[i][j][k][14];	  
	  if (k==subev) h_maqe_sub->SetBinContent(e+1, i+1, maqe_sub[i][j][k][e]);
	}
	g_maqe_sub[i][j][k] = new TSpline3(Form("maqe_p%d_e%d_ev%d_sub",i,j,k+1), maqe, maqe_sub[i][j][k], np_maqe);
	g_pfo_sub[i][j][k] = new TSpline3(Form("pfo_p%d_e%d_ev%d_sub",i,j,k+1), pfo, pfo_sub[i][j][k], np);
	g_ebo_sub[i][j][k] = new TSpline3(Form("ebo_p%d_e%d_ev%d_sub",i,j,k+1), ebo, ebo_sub[i][j][k], np);
	
	if (cc1pi_sub[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    ca5_cc1pi_sub[i][j][k][e] /= cc1pi_sub[i][j][k];
	    manff_cc1pi_sub[i][j][k][e] /= cc1pi_sub[i][j][k];
	    bgscl_cc1pi_sub[i][j][k][e] /= cc1pi_sub[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    ca5_cc1pi_sub[i][j][k][e] = 1;
	    manff_cc1pi_sub[i][j][k][e] = 1;
	    bgscl_cc1pi_sub[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  ca5_cc1pi_sub[i][j][k][e] /= ca5_cc1pi_sub[i][j][k][6];
	  manff_cc1pi_sub[i][j][k][e] /= manff_cc1pi_sub[i][j][k][6];
	  bgscl_cc1pi_sub[i][j][k][e] /= bgscl_cc1pi_sub[i][j][k][6];
	  if (k==subev) {
	    h_ca5_cc1pi_sub->SetBinContent(e+1, i+1, ca5_cc1pi_sub[i][j][k][e]);
	    h_manff_cc1pi_sub->SetBinContent(e+1, i+1, manff_cc1pi_sub[i][j][k][e]);
	    h_bgscl_cc1pi_sub->SetBinContent(e+1, i+1, bgscl_cc1pi_sub[i][j][k][e]);
	  }
	}
	g_ca5_cc1pi_sub[i][j][k] = new TSpline3(Form("ca5_cc1pi_p%d_e%d_ev%d_sub",i,j,k+1), ca5_cc1pi, ca5_cc1pi_sub[i][j][k], np);
	g_manff_cc1pi_sub[i][j][k] = new TSpline3(Form("manff_cc1pi_p%d_e%d_ev%d_sub",i,j,k+1), manff_cc1pi, manff_cc1pi_sub[i][j][k], np);
	g_bgscl_cc1pi_sub[i][j][k] = new TSpline3(Form("bgscl_cc1pi_p%d_e%d_ev%d_sub",i,j,k+1), bgscl_cc1pi, bgscl_cc1pi_sub[i][j][k], np);
	
	if (ncpiz_sub[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpiz_sub[i][j][k][e] /= ncpiz_sub[i][j][k];
	    manff_ncpiz_sub[i][j][k][e] /= ncpiz_sub[i][j][k];
	    bgscl_ncpiz_sub[i][j][k][e] /= ncpiz_sub[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpiz_sub[i][j][k][e] = 1;
	    manff_ncpiz_sub[i][j][k][e] = 1;
	    bgscl_ncpiz_sub[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  ca5_ncpiz_sub[i][j][k][e] /= ca5_ncpiz_sub[i][j][k][6];
	  manff_ncpiz_sub[i][j][k][e] /= manff_ncpiz_sub[i][j][k][6];
	  bgscl_ncpiz_sub[i][j][k][e] /= bgscl_ncpiz_sub[i][j][k][6];
	  if (k==subev) {
	    h_ca5_ncpiz_sub->SetBinContent(e+1, i+1, ca5_ncpiz_sub[i][j][k][e]);
	    h_manff_ncpiz_sub->SetBinContent(e+1, i+1, manff_ncpiz_sub[i][j][k][e]);
	    h_bgscl_ncpiz_sub->SetBinContent(e+1, i+1, bgscl_ncpiz_sub[i][j][k][e]);
	  }
	}
	g_ca5_ncpiz_sub[i][j][k] = new TSpline3(Form("ca5_ncpiz_p%d_e%d_ev%d_sub",i,j,k+1), ca5_ncpiz, ca5_ncpiz_sub[i][j][k], np);
	g_manff_ncpiz_sub[i][j][k] = new TSpline3(Form("manff_ncpiz_p%d_e%d_ev%d_sub",i,j,k+1), manff_ncpiz, manff_ncpiz_sub[i][j][k], np);
	g_bgscl_ncpiz_sub[i][j][k] = new TSpline3(Form("bgscl_ncpiz_p%d_e%d_ev%d_sub",i,j,k+1), bgscl_ncpiz, bgscl_ncpiz_sub[i][j][k], np);

	if (ncpipm_sub[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpipm_sub[i][j][k][e] /= ncpipm_sub[i][j][k];
	    manff_ncpipm_sub[i][j][k][e] /= ncpipm_sub[i][j][k];
	    bgscl_ncpipm_sub[i][j][k][e] /= ncpipm_sub[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpipm_sub[i][j][k][e] = 1;
	    manff_ncpipm_sub[i][j][k][e] = 1;
	    bgscl_ncpipm_sub[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  ca5_ncpipm_sub[i][j][k][e] /= ca5_ncpipm_sub[i][j][k][6];
	  manff_ncpipm_sub[i][j][k][e] /= manff_ncpipm_sub[i][j][k][6];
	  bgscl_ncpipm_sub[i][j][k][e] /= bgscl_ncpipm_sub[i][j][k][6];
	  if (k==subev) {
	    h_ca5_ncpipm_sub->SetBinContent(e+1, i+1, ca5_ncpipm_sub[i][j][k][e]);
	    h_manff_ncpipm_sub->SetBinContent(e+1, i+1, manff_ncpipm_sub[i][j][k][e]);
	    h_bgscl_ncpipm_sub->SetBinContent(e+1, i+1, bgscl_ncpipm_sub[i][j][k][e]);
	  }
	}
	g_ca5_ncpipm_sub[i][j][k] = new TSpline3(Form("ca5_ncpipm_p%d_e%d_ev%d_sub",i,j,k+1), ca5_ncpipm, ca5_ncpipm_sub[i][j][k], np);
	g_manff_ncpipm_sub[i][j][k] = new TSpline3(Form("manff_ncpipm_p%d_e%d_ev%d_sub",i,j,k+1), manff_ncpipm, manff_ncpipm_sub[i][j][k], np);
	g_bgscl_ncpipm_sub[i][j][k] = new TSpline3(Form("bgscl_ncpipm_p%d_e%d_ev%d_sub",i,j,k+1), bgscl_ncpipm, bgscl_ncpipm_sub[i][j][k], np);       

	if (ccoth_sub[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    dismpishp_sub[i][j][k][e] /= ccoth_sub[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    dismpishp_sub[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  dismpishp_sub[i][j][k][e] /= dismpishp_sub[i][j][k][6];
	  if (k==subev) h_dismpishp_sub->SetBinContent(e+1, i+1, dismpishp_sub[i][j][k][e]);
	}
	g_dismpishp_sub[i][j][k] = new TSpline3(Form("dismpishp_p%d_e%d_ev%d_sub",i,j,k+1), dismpishp, dismpishp_sub[i][j][k], np);

	g_maqe_sub[i][j][k]->SetName(Form("maqe_p%d_e%d_ev%d_sub",i,j,k+1));
	g_pfo_sub[i][j][k]->SetName(Form("pfo_p%d_e%d_ev%d_sub",i,j,k+1));
	g_ebo_sub[i][j][k]->SetName(Form("ebo_p%d_e%d_ev%d_sub",i,j,k+1));
	g_ca5_cc1pi_sub[i][j][k]->SetName(Form("ca5_cc1pi_p%d_e%d_ev%d_sub",i,j,k+1));
	g_ca5_ncpiz_sub[i][j][k]->SetName(Form("ca5_ncpiz_p%d_e%d_ev%d_sub",i,j,k+1));
	g_ca5_ncpipm_sub[i][j][k]->SetName(Form("ca5_ncpipm_p%d_e%d_ev%d_sub",i,j,k+1));
	g_manff_cc1pi_sub[i][j][k]->SetName(Form("manff_cc1pi_p%d_e%d_ev%d_sub",i,j,k+1));
	g_manff_ncpiz_sub[i][j][k]->SetName(Form("manff_ncpiz_p%d_e%d_ev%d_sub",i,j,k+1));
	g_manff_ncpipm_sub[i][j][k]->SetName(Form("manff_ncpipm_p%d_e%d_ev%d_sub",i,j,k+1));
	g_bgscl_cc1pi_sub[i][j][k]->SetName(Form("bgscl_cc1pi_p%d_e%d_ev%d_sub",i,j,k+1));
	g_bgscl_ncpiz_sub[i][j][k]->SetName(Form("bgscl_ncpiz_p%d_e%d_ev%d_sub",i,j,k+1));
	g_bgscl_ncpipm_sub[i][j][k]->SetName(Form("bgscl_ncpipm_p%d_e%d_ev%d_sub",i,j,k+1));
	g_dismpishp_sub[i][j][k]->SetName(Form("dismpishp_p%d_e%d_ev%d_sub",i,j,k+1));

	outfile->cd();
	g_maqe_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_pfo_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ebo_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ca5_cc1pi_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ca5_ncpiz_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ca5_ncpipm_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_manff_cc1pi_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_manff_ncpiz_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_manff_ncpipm_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_bgscl_cc1pi_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_bgscl_ncpiz_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_bgscl_ncpipm_sub[i][j][k]->Write(0, TObject::kWriteDelete);
	g_dismpishp_sub[i][j][k]->Write(0, TObject::kWriteDelete);	
      }
    }
  }

  // construct and save splines to file -- multi-GeV
  for (int i = 0; i < pidmultbins; ++i) {
    for (int j = 1; j < 2/*etruemultbins*/; ++j) {
      for (int k = 0; k < 2; ++k) {

	if (ccqe_mult[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    maqe_mult[i][j][k][e] /= ccqe_mult[i][j][k];
	    pfo_mult[i][j][k][e] /= ccqe_mult[i][j][k];
	    ebo_mult[i][j][k][e] /= ccqe_mult[i][j][k];
	  }
	  for (int e = np; e < np_maqe; ++e) {
	    maqe_mult[i][j][k][e] /= ccqe_mult[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    maqe_mult[i][j][k][e] = 1;
	    pfo_mult[i][j][k][e] = 1;
	    ebo_mult[i][j][k][e] = 1;
	  }
	  for (int e = np; e < np_maqe; ++e) {
	    maqe_mult[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  maqe_mult[i][j][k][e] /= maqe_mult[i][j][k][14];
	  pfo_mult[i][j][k][e] /= pfo_mult[i][j][k][6];
	  ebo_mult[i][j][k][e] /= ebo_mult[i][j][k][6];
	  if (k==subev) {
	    h_maqe_mult->SetBinContent(e+1, i+1, maqe_mult[i][j][k][e]);
	    h_pfo_mult->SetBinContent(e+1, i+1, pfo_mult[i][j][k][e]);
	    h_ebo_mult->SetBinContent(e+1, i+1, ebo_mult[i][j][k][e]);
	  }
	}
	for (int e = np; e < np_maqe; ++e) {
	  maqe_mult[i][j][k][e] /= maqe_mult[i][j][k][14];
	  if (k==subev) h_maqe_mult->SetBinContent(e+1, i+1, maqe_mult[i][j][k][e]);
	}
	g_maqe_mult[i][j][k] = new TSpline3(Form("maqe_p%d_e%d_ev%d_mult",i,j,k+1), maqe, maqe_mult[i][j][k], np_maqe);
	g_pfo_mult[i][j][k] = new TSpline3(Form("pfo_p%d_e%d_ev%d_mult",i,j,k+1), pfo, pfo_mult[i][j][k], np);
	g_ebo_mult[i][j][k] = new TSpline3(Form("ebo_p%d_e%d_ev%d_mult",i,j,k+1), ebo, ebo_mult[i][j][k], np);

	if (cc1pi_mult[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    ca5_cc1pi_mult[i][j][k][e] /= cc1pi_mult[i][j][k];
	    manff_cc1pi_mult[i][j][k][e] /= cc1pi_mult[i][j][k];
	    bgscl_cc1pi_mult[i][j][k][e] /= cc1pi_mult[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    ca5_cc1pi_mult[i][j][k][e] = 1;
	    manff_cc1pi_mult[i][j][k][e] = 1;
	    bgscl_cc1pi_mult[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  ca5_cc1pi_mult[i][j][k][e] /= ca5_cc1pi_mult[i][j][k][6];
	  manff_cc1pi_mult[i][j][k][e] /= manff_cc1pi_mult[i][j][k][6];
	  bgscl_cc1pi_mult[i][j][k][e] /= bgscl_cc1pi_mult[i][j][k][6];
	  if (k==subev) {
	    h_ca5_cc1pi_mult->SetBinContent(e+1, i+1, ca5_cc1pi_mult[i][j][k][e]);
	    h_manff_cc1pi_mult->SetBinContent(e+1, i+1, manff_cc1pi_mult[i][j][k][e]);
	    h_bgscl_cc1pi_mult->SetBinContent(e+1, i+1, bgscl_cc1pi_mult[i][j][k][e]);
	  }
	}
	g_ca5_cc1pi_mult[i][j][k] = new TSpline3(Form("ca5_cc1pi_p%d_e%d_ev%d_mult",i,j,k+1), ca5_cc1pi, ca5_cc1pi_mult[i][j][k], np);
	g_manff_cc1pi_mult[i][j][k] = new TSpline3(Form("manff_cc1pi_p%d_e%d_ev%d_mult",i,j,k+1), manff_cc1pi, manff_cc1pi_mult[i][j][k], np);
	g_bgscl_cc1pi_mult[i][j][k] = new TSpline3(Form("bgscl_cc1pi_p%d_e%d_ev%d_mult",i,j,k+1), bgscl_cc1pi, bgscl_cc1pi_mult[i][j][k], np);
	
	if (ncpiz_mult[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpiz_mult[i][j][k][e] /= ncpiz_mult[i][j][k];
	    manff_ncpiz_mult[i][j][k][e] /= ncpiz_mult[i][j][k];
	    bgscl_ncpiz_mult[i][j][k][e] /= ncpiz_mult[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpiz_mult[i][j][k][e] = 1;
	    manff_ncpiz_mult[i][j][k][e] = 1;
	    bgscl_ncpiz_mult[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  ca5_ncpiz_mult[i][j][k][e] /= ca5_ncpiz_mult[i][j][k][6];
	  manff_ncpiz_mult[i][j][k][e] /= manff_ncpiz_mult[i][j][k][6];
	  bgscl_ncpiz_mult[i][j][k][e] /= bgscl_ncpiz_mult[i][j][k][6];
	  if (k==subev) {
	    h_ca5_ncpiz_mult->SetBinContent(e+1, i+1, ca5_ncpiz_mult[i][j][k][e]);
	    h_manff_ncpiz_mult->SetBinContent(e+1, i+1, manff_ncpiz_mult[i][j][k][e]);
	    h_bgscl_ncpiz_mult->SetBinContent(e+1, i+1, bgscl_ncpiz_mult[i][j][k][e]);
	  }
	}
	g_ca5_ncpiz_mult[i][j][k] = new TSpline3(Form("ca5_ncpiz_p%d_e%d_ev%d_mult",i,j,k+1), ca5_ncpiz, ca5_ncpiz_mult[i][j][k], np);
	g_manff_ncpiz_mult[i][j][k] = new TSpline3(Form("manff_ncpiz_p%d_e%d_ev%d_mult",i,j,k+1), manff_ncpiz, manff_ncpiz_mult[i][j][k], np);
	g_bgscl_ncpiz_mult[i][j][k] = new TSpline3(Form("bgscl_ncpiz_p%d_e%d_ev%d_mult",i,j,k+1), bgscl_ncpiz, bgscl_ncpiz_mult[i][j][k], np);

	if (ncpipm_mult[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpipm_mult[i][j][k][e] /= ncpipm_mult[i][j][k];
	    manff_ncpipm_mult[i][j][k][e] /= ncpipm_mult[i][j][k];
	    bgscl_ncpipm_mult[i][j][k][e] /= ncpipm_mult[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    ca5_ncpipm_mult[i][j][k][e] = 1;
	    manff_ncpipm_mult[i][j][k][e] = 1;
	    bgscl_ncpipm_mult[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  ca5_ncpipm_mult[i][j][k][e] /= ca5_ncpipm_mult[i][j][k][6];
	  manff_ncpipm_mult[i][j][k][e] /= manff_ncpipm_mult[i][j][k][6];
	  bgscl_ncpipm_mult[i][j][k][e] /= bgscl_ncpipm_mult[i][j][k][6];
	  if (k==subev) {
	    h_ca5_ncpipm_mult->SetBinContent(e+1, i+1, ca5_ncpipm_mult[i][j][k][e]);
	    h_manff_ncpipm_mult->SetBinContent(e+1, i+1, manff_ncpipm_mult[i][j][k][e]);
	    h_bgscl_ncpipm_mult->SetBinContent(e+1, i+1, bgscl_ncpipm_mult[i][j][k][e]);
	  }
	}
	g_ca5_ncpipm_mult[i][j][k] = new TSpline3(Form("ca5_ncpipm_p%d_e%d_ev%d_mult",i,j,k+1), ca5_ncpipm, ca5_ncpipm_mult[i][j][k], np);
	g_manff_ncpipm_mult[i][j][k] = new TSpline3(Form("manff_ncpipm_p%d_e%d_ev%d_mult",i,j,k+1), manff_ncpipm, manff_ncpipm_mult[i][j][k], np);
	g_bgscl_ncpipm_mult[i][j][k] = new TSpline3(Form("bgscl_ncpipm_p%d_e%d_ev%d_mult",i,j,k+1), bgscl_ncpipm, bgscl_ncpipm_mult[i][j][k], np);       

	if (ccoth_mult[i][j][k]>1e-3) {
	  for (int e = 0; e < np; ++e) {
	    dismpishp_mult[i][j][k][e] /= ccoth_mult[i][j][k];
	  }
	} else {
	  for (int e = 0; e < np; ++e) {
	    dismpishp_mult[i][j][k][e] = 1;
	  }
	}
	for (int e = 0; e < np; ++e) {
	  dismpishp_mult[i][j][k][e] /= dismpishp_mult[i][j][k][6];
	  if (k==subev) h_dismpishp_mult->SetBinContent(e+1, i+1, dismpishp_mult[i][j][k][e]);
	}
	g_dismpishp_mult[i][j][k] = new TSpline3(Form("dismpishp_p%d_e%d_ev%d_mult",i,j,k+1), dismpishp, dismpishp_mult[i][j][k], np);

	g_maqe_mult[i][j][k]->SetName(Form("maqe_p%d_e%d_ev%d_mult",i,j,k+1));
	g_pfo_mult[i][j][k]->SetName(Form("pfo_p%d_e%d_ev%d_mult",i,j,k+1));
	g_ebo_mult[i][j][k]->SetName(Form("ebo_p%d_e%d_ev%d_mult",i,j,k+1));
	g_ca5_cc1pi_mult[i][j][k]->SetName(Form("ca5_cc1pi_p%d_e%d_ev%d_mult",i,j,k+1));
	g_ca5_ncpiz_mult[i][j][k]->SetName(Form("ca5_ncpiz_p%d_e%d_ev%d_mult",i,j,k+1));
	g_ca5_ncpipm_mult[i][j][k]->SetName(Form("ca5_ncpipm_p%d_e%d_ev%d_mult",i,j,k+1));
	g_manff_cc1pi_mult[i][j][k]->SetName(Form("manff_cc1pi_p%d_e%d_ev%d_mult",i,j,k+1));
	g_manff_ncpiz_mult[i][j][k]->SetName(Form("manff_ncpiz_p%d_e%d_ev%d_mult",i,j,k+1));
	g_manff_ncpipm_mult[i][j][k]->SetName(Form("manff_ncpipm_p%d_e%d_ev%d_mult",i,j,k+1));
	g_bgscl_cc1pi_mult[i][j][k]->SetName(Form("bgscl_cc1pi_p%d_e%d_ev%d_mult",i,j,k+1));
	g_bgscl_ncpiz_mult[i][j][k]->SetName(Form("bgscl_ncpiz_p%d_e%d_ev%d_mult",i,j,k+1));
	g_bgscl_ncpipm_mult[i][j][k]->SetName(Form("bgscl_ncpipm_p%d_e%d_ev%d_mult",i,j,k+1));
	g_dismpishp_mult[i][j][k]->SetName(Form("dismpishp_p%d_e%d_ev%d_mult",i,j,k+1));

	outfile->cd();
	g_maqe_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_pfo_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ebo_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ca5_cc1pi_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ca5_ncpiz_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_ca5_ncpipm_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_manff_cc1pi_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_manff_ncpiz_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_manff_ncpipm_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_bgscl_cc1pi_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_bgscl_ncpiz_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_bgscl_ncpipm_mult[i][j][k]->Write(0, TObject::kWriteDelete);
	g_dismpishp_mult[i][j][k]->Write(0, TObject::kWriteDelete);	
      }
    }
  }

  outfile->Close();

  plotfile->cd();
  h_maqe_sub->Write(0, TObject::kWriteDelete);// one sub-event and two sub-events
  h_maqe_mult->Write(0, TObject::kWriteDelete);
  h_pfo_sub->Write(0, TObject::kWriteDelete);
  h_pfo_mult->Write(0, TObject::kWriteDelete);
  h_ebo_sub->Write(0, TObject::kWriteDelete);
  h_ebo_mult->Write(0, TObject::kWriteDelete);
  h_ca5_cc1pi_sub->Write(0, TObject::kWriteDelete);
  h_ca5_cc1pi_mult->Write(0, TObject::kWriteDelete);
  h_ca5_ncpiz_sub->Write(0, TObject::kWriteDelete);
  h_ca5_ncpiz_mult->Write(0, TObject::kWriteDelete);
  h_ca5_ncpipm_sub->Write(0, TObject::kWriteDelete);
  h_ca5_ncpipm_mult->Write(0, TObject::kWriteDelete);
  h_manff_cc1pi_sub->Write(0, TObject::kWriteDelete);
  h_manff_cc1pi_mult->Write(0, TObject::kWriteDelete);
  h_manff_ncpiz_sub->Write(0, TObject::kWriteDelete);
  h_manff_ncpiz_mult->Write(0, TObject::kWriteDelete);
  h_manff_ncpipm_sub->Write(0, TObject::kWriteDelete);
  h_manff_ncpipm_mult->Write(0, TObject::kWriteDelete);
  h_bgscl_cc1pi_sub->Write(0, TObject::kWriteDelete);
  h_bgscl_cc1pi_mult->Write(0, TObject::kWriteDelete);
  h_bgscl_ncpiz_sub->Write(0, TObject::kWriteDelete);
  h_bgscl_ncpiz_mult->Write(0, TObject::kWriteDelete);
  h_bgscl_ncpipm_sub->Write(0, TObject::kWriteDelete);
  h_bgscl_ncpipm_mult->Write(0, TObject::kWriteDelete);
  h_dismpishp_sub->Write(0, TObject::kWriteDelete);
  h_dismpishp_mult->Write(0, TObject::kWriteDelete);
  plotfile->Close();
  
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
