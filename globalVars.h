#include <string>
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TBranch.h"

  static const std::string filename = "jan14sk4_skdetsim13p90_neut532.reduc.fQv4r0.root";
  static const std::string splfilename = "/storage/shared/xiaoyueli/skmc/test_fc/splines/jan14sk4_skdetsim13p90_neut532.reduc.fQv4r.splines.root";
  static const std::string xsecfilename = "xsec_covariance_2015a_q3_1.2_withRPA_v1.root";
  static const std::string splinepath = "/storage/shared/xiaoyueli/skmc/test_fc/binned_splines";
  static const std::string evsplnames[] = {
    "byEv_maqe_ccqe_gr",
    "byEv_pfo_ccqe_gr",
    "byEv_ebo_ccqe_gr",
    "byEv_ca5_cc1pi_gr",
    "byEv_ca5_ncpiz_gr",
    "byEv_ca5_ncpipm_gr",
    "byEv_manff_cc1pi_gr",
    "byEv_manff_ncpiz_g ",
    "byEv_manff_ncpipm_gr",
    "byEv_bgscl_cc1pi_gr",
    "byEv_bgscl_ncpiz_gr",
    "byEv_bgscl_ncpipm_gr",
    "byEv_dismpishp_ccoth_gr",
    "byEv_rpa_ccqe_gr",
    "byEv_sccvec_ccqe_gr",
    "byEv_sccvec_ncoth_gr",
    "byEv_sccaxl_ccqe_gr",
    "byEv_sccaxl_ncoth_gr"
};

  static const std::string intnames[] = {
    "ccqe", // mode = 1; 0
    "cc1pi", // mode = [11,13]; 1
    "cccoh", // mode = 16; 2
    "ccoth", // mode = [17,30); 3
    "ncpiz", // mode = [31,32]; 4
    "ncpipm", // mode = [33,34]; 5
    "nccoh", // mode = 36; 6
    "ncoth", // mode = [37,52]; 7
    "mec" // mode = 2; 8
  };

  static const double scale = 1993.6 / 99. /365.; // live-time 1993.6 days

  Int_t nring, fqnse, nmue, evclass, ipnu[50], mode, fqmrnring[20];
  UInt_t ip[50], nhitac;
  Float_t evis, dir[10][3], wall, amomm[10], amome[10], wallv, pnu[50], fq1rmom[10][7], fq1rnll[10][7], fq1rdir[10][7][3], oscwgt;
  Double_t wgtosc, wgtflx, fqdwall, coszenith, coslep;


  TGraph          *byEv_maqe_ccqe_gr;
  TGraph          *byEv_pfo_ccqe_gr;
  TGraph          *byEv_ebo_ccqe_gr;
  TGraph          *byEv_ca5_cc1pi_gr;
  TGraph          *byEv_ca5_ncpiz_gr;
  TGraph          *byEv_ca5_ncpipm_gr;
  TGraph          *byEv_manff_cc1pi_gr;
  TGraph          *byEv_manff_ncpiz_gr;
  TGraph          *byEv_manff_ncpipm_gr;
  TGraph          *byEv_bgscl_cc1pi_gr;
  TGraph          *byEv_bgscl_ncpiz_gr;
  TGraph          *byEv_bgscl_ncpipm_gr;
  TGraph          *byEv_dismpishp_ccoth_gr;
  TGraph          *byEv_rpa_ccqe_gr;
  //TGraph          *byEv_sccvec_ccqe_gr;
  //TGraph          *byEv_sccvec_ncoth_gr;
  //TGraph          *byEv_sccaxl_ccqe_gr;
  //TGraph          *byEv_sccaxl_ncoth_gr;
  TBranch          *byEv_maqe_ccqe_br;
  TBranch          *byEv_pfo_ccqe_br;
  TBranch          *byEv_ebo_ccqe_br;
  TBranch          *byEv_ca5_cc1pi_br;
  TBranch          *byEv_ca5_ncpiz_br;
  TBranch          *byEv_ca5_ncpipm_br;
  TBranch          *byEv_manff_cc1pi_br;
  TBranch          *byEv_manff_ncpiz_br;
  TBranch          *byEv_manff_ncpipm_br;
  TBranch          *byEv_bgscl_cc1pi_br;
  TBranch          *byEv_bgscl_ncpiz_br;
  TBranch          *byEv_bgscl_ncpipm_br;
  TBranch          *byEv_dismpishp_ccoth_br;
  TBranch          *byEv_rpa_ccqe_br;
  //TBranch          *byEv_sccvec_ccqe_br;
  //TBranch          *byEv_sccvec_ncoth_br;
  //TBranch          *byEv_sccaxl_ccqe_br;
  //TBranch          *byEv_sccaxl_ncoth_br;

  TFile *nuespline;
  TFile *nuebarspline;
  TFile *numuspline;
  TFile *numubarspline;

  TH1D *dmaqe[4][2];
  TH1D *dpfo[4][2];
  TH1D *debo[4][2];
  TH1D *dca5_cc1pi[4][2];
  TH1D *dca5_ncpiz[4][2];
  TH1D *dca5_ncpipm[4][2];
  TH1D *dmanff_cc1pi[4][2];
  TH1D *dmanff_ncpiz[4][2];
  TH1D *dmanff_ncpipm[4][2];
  TH1D *dbgscl_cc1pi[4][2];
  TH1D *dbgscl_ncpiz[4][2];
  TH1D *dbgscl_ncpipm[4][2];
  TH1D *ddismpishp[4][2];
