#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TArrayF.h"
#include "TGraph.h"
#include "TChain.h"

std::string fInputFileName;
std::string fOutputFileName;
std::string fInputSKFileName;
bool fUseWildcard;

static const std::string branchnames[13] = {
  "byEv_maqe_ccqe_gr",
  "byEv_pfo_ccqe_gr",
  "byEv_ebo_ccqe_gr",
  "byEv_ca5_cc1pi_gr",
  "byEv_ca5_ncpiz_gr",
  "byEv_ca5_ncpipm_gr",
  "byEv_manff_cc1pi_gr",
  "byEv_manff_ncpiz_gr",
  "byEv_manff_ncpipm_gr",
  "byEv_bgscl_cc1pi_gr",
  "byEv_bgscl_ncpiz_gr",
  "byEv_bgscl_ncpipm_gr",
  "byEv_dismpishp_ccoth_gr"
};

void Usage();
int ParseArgs(int argc, char *argv[]);

int main(int argc, char * argv[])
{
  int args = ParseArgs(argc, argv);
  if(args) {
    std::cerr << "Usage: " << std::endl;
    Usage();
  }
  TChain *sktree;
  TChain *wtree;
  if(!fUseWildcard) {
    TFile *infile = new TFile(fInputFileName.c_str(), "read");
    if(!infile) {
      std::cerr << "Cannot open T2KReWeight-style input file!" << std::endl;
      exit(1);
    }
    wtree = (TChain*)infile->Get("weightstree");
    if(!wtree) {
      std::cerr << "Cannot find weightstree!" << std::endl;
      exit(1);
    }
    TFile *skfile = new TFile(fInputSKFileName.c_str(), "read");
    if(!skfile) {
      std::cerr << "Cannot open SK input file!" << std::endl;
      exit(1);
    }
    sktree = (TChain*)skfile->Get("h1");
    if(!sktree) {
      std::cerr << "Cannot find h1!" << std::endl;
      exit(1);
    }
  } else {
    sktree = new TChain("h1");
    int nskfiles = sktree->Add(fInputSKFileName.c_str());
    if(!nskfiles) {
      std::cerr << "Cannot find SK input files!" <<std::endl;
      exit(1);
    }
    wtree = new TChain("weightstree");
    int nwgtfiles = wtree->Add(fInputFileName.c_str());
    if(!nwgtfiles) {
      std::cerr << "Cannot find T2KReWeight input files!" << std::endl;
      exit(1);
    }
  }

  int mode, nweights;
  TArrayF *weights = 0;
  TArrayF *brNXSec_MaCCQE = 0;
  TArrayF *brNIWG2014a_pF_O16 = 0;
  TArrayF *brNIWG2014a_Eb_O16 = 0;
  TArrayF *brNXSec_CA5RES = 0;
  TArrayF *brNXSec_MaNFFRES = 0;
  TArrayF *brNXSec_BgSclRES = 0;
  TArrayF *brNIWG2012a_dismpishp = 0;

  sktree->SetBranchAddress("mode", &mode);
  wtree->SetBranchAddress("nweights", &nweights);
  wtree->SetBranchAddress("weights", &weights);
  wtree->SetBranchAddress("brNXSec_MaCCQE", &brNXSec_MaCCQE);
  wtree->SetBranchAddress("brNIWG2014a_pF_O16", &brNIWG2014a_pF_O16);
  wtree->SetBranchAddress("brNIWG2014a_Eb_O16", &brNIWG2014a_Eb_O16);
  wtree->SetBranchAddress("brNXSec_CA5RES", &brNXSec_CA5RES);
  wtree->SetBranchAddress("brNXSec_MaNFFRES", &brNXSec_MaNFFRES);
  wtree->SetBranchAddress("brNXSec_BgSclRES", &brNXSec_BgSclRES);
  wtree->SetBranchAddress("brNIWG2012a_dismpishp", &brNIWG2012a_dismpishp);
  sktree->AddFriend(wtree);
  double maqe[23], pfo[13], ebo[13], ca5[13], manff[13], bgscl[13], dismpishp[13];
  int maqenom = 14;
  int pfonom = 21 + 6;
  int ebonom = 21 + 13 + 6;
  int ca5nom = 21 + 13 + 13 + 6;
  int manffnom = 21 + 13 + 13 + 13 + 6;
  int bgsclnom = 21 + 13 + 13 + 13 + 13 + 6;
  int dismpishpnom = 21 + 13 + 13 + 13 + 13 + 13 + 6;

  wtree->GetEntry(0);
  for(int i=0; i<13; ++i) {
    maqe[i] = brNXSec_MaCCQE->At(i);
    pfo[i] = brNIWG2014a_pF_O16->At(21+i);
    ebo[i] = brNIWG2014a_Eb_O16->At(21+13+i);
    ca5[i] = brNXSec_CA5RES->At(21+13+13+i);
    manff[i] = brNXSec_MaNFFRES->At(21+13+13+13+i);
    bgscl[i] = brNXSec_BgSclRES->At(21+13+13+13+13+i);
    dismpishp[i] = brNIWG2012a_dismpishp->At(21+13+13+13+13+13+i);
  }
  for(int i=13; i<21; ++i) maqe[i] = brNXSec_MaCCQE->At(i);

  TFile *outfile = new TFile(fOutputFileName.c_str(), "recreate");
  TTree *spltree = new TTree("SplinesByEvent", "event-by-event spline");
  TGraph *byEv_maqe_ccqe_gr = new TGraph(21);
  byEv_maqe_ccqe_gr->SetName("byEv_maqe_ccqe_gr");
  TGraph *byEv_pfo_ccqe_gr = new TGraph(13);
  byEv_pfo_ccqe_gr->SetName("byEv_pfo_ccqe_gr");
  TGraph *byEv_ebo_ccqe_gr = new TGraph(13);
  byEv_ebo_ccqe_gr->SetName("byEv_ebo_ccqe_gr");
  TGraph *byEv_ca5_cc1pi_gr = new TGraph(13);
  byEv_ca5_cc1pi_gr->SetName("byEv_ca5_cc1pi_gr");
  TGraph *byEv_ca5_ncpiz_gr = new TGraph(13);
  byEv_ca5_ncpiz_gr->SetName("byEv_ca5_ncpiz_gr");
  TGraph *byEv_ca5_ncpipm_gr = new TGraph(13);
  byEv_ca5_ncpipm_gr->SetName("byEv_ca5_ncpipm_gr");
  TGraph *byEv_manff_cc1pi_gr = new TGraph(13);
  byEv_manff_cc1pi_gr->SetName("byEv_manff_cc1pi_gr");
  TGraph *byEv_manff_ncpiz_gr = new TGraph(13);
  byEv_manff_ncpiz_gr->SetName("byEv_manff_ncpiz_gr");
  TGraph *byEv_manff_ncpipm_gr = new TGraph(13);
  byEv_manff_ncpipm_gr->SetName("byEv_manff_ncpipm_gr");
  TGraph *byEv_bgscl_cc1pi_gr = new TGraph(13);
  byEv_bgscl_cc1pi_gr->SetName("byEv_bgscl_cc1pi_gr");
  TGraph *byEv_bgscl_ncpiz_gr = new TGraph(13);
  byEv_bgscl_ncpiz_gr->SetName("byEv_bgscl_ncpiz_gr");
  TGraph *byEv_bgscl_ncpipm_gr = new TGraph(13);
  byEv_bgscl_ncpipm_gr->SetName("byEv_bgscl_ncpipm_gr");
  TGraph *byEv_dismpishp_ccoth_gr = new TGraph(13);
  byEv_dismpishp_ccoth_gr->SetName("byEv_dismpishp_ccoth_gr");
  spltree->Branch("byEv_maqe_ccqe_gr", "TGraph", &byEv_maqe_ccqe_gr, 1280000, 1);
  spltree->Branch("byEv_pfo_ccqe_gr", "TGraph", &byEv_pfo_ccqe_gr, 1280000, 1);
  spltree->Branch("byEv_ebo_ccqe_gr", "TGraph", &byEv_ebo_ccqe_gr, 1280000, 1);
  spltree->Branch("byEv_ca5_cc1pi_gr", "TGraph", &byEv_ca5_cc1pi_gr, 1280000, 1);
  spltree->Branch("byEv_ca5_ncpiz_gr", "TGraph", &byEv_ca5_ncpiz_gr, 1280000, 1);
  spltree->Branch("byEv_ca5_ncpipm_gr", "TGraph", &byEv_ca5_ncpipm_gr, 1280000, 1);
  spltree->Branch("byEv_manff_cc1pi_gr", "TGraph", &byEv_manff_cc1pi_gr, 1280000, 1);
  spltree->Branch("byEv_manff_ncpiz_gr", "TGraph", &byEv_manff_ncpiz_gr, 1280000, 1);
  spltree->Branch("byEv_manff_ncpipm_gr", "TGraph", &byEv_manff_ncpipm_gr, 1280000, 1);
  spltree->Branch("byEv_bgscl_cc1pi_gr", "TGraph", &byEv_bgscl_cc1pi_gr, 1280000, 1);
  spltree->Branch("byEv_bgscl_ncpiz_gr", "TGraph", &byEv_bgscl_ncpiz_gr, 1280000, 1);
  spltree->Branch("byEv_bgscl_ncpipm_gr", "TGraph", &byEv_bgscl_ncpipm_gr, 1280000, 1);
  spltree->Branch("byEv_dismpishp_ccoth_gr", "TGraph", &byEv_dismpishp_ccoth_gr, 1280000, 1);

  int nentries = sktree->GetEntries();
  if(wtree->GetEntries() != nentries) {
    std::cerr << "sktree and weightstree don't have the same length!" << std::endl;
    exit(1);
  }
  for(int i = 0; i < nentries; ++i) {
    sktree->GetEntry(i);
    wtree->GetEntry(i);

    if(abs(mode)==1) { // CCQE: maqe, pfo, ebo
      for(int j=0; j<13; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], (weights->At(maqenom) > 0 ? weights->At(j)/weights->At(maqenom) : 1.));
	byEv_pfo_ccqe_gr->SetPoint(j, pfo[j], (weights->At(pfonom) > 0 ? weights->At(21+j)/weights->At(pfonom) : 1.));
	byEv_ebo_ccqe_gr->SetPoint(j, ebo[j], (weights->At(ebonom) > 0 ? weights->At(34+j)/weights->At(ebonom) : 1.));
	byEv_ca5_cc1pi_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpiz_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpipm_gr->SetPoint(j, ca5[j], 1);
	byEv_manff_cc1pi_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpiz_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpipm_gr->SetPoint(j, manff[j], 1);
	byEv_bgscl_cc1pi_gr->SetPoint(j, bgscl[j], 1);
	byEv_bgscl_ncpiz_gr->SetPoint(j, bgscl[j], 1);
	byEv_bgscl_ncpipm_gr->SetPoint(j, bgscl[j], 1);
	byEv_dismpishp_ccoth_gr->SetPoint(j, dismpishp[j], 1);
      }
      for(int j=13; j<21; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], (weights->At(maqenom) > 0 ? weights->At(j)/weights->At(maqenom) : 1.));
      }
    }

    else if(abs(mode) >=11 && abs(mode) <=13) { // cc1pi: ca5, manff, bgscl
      for(int j=0; j<13; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
	byEv_pfo_ccqe_gr->SetPoint(j, pfo[j], 1);
	byEv_ebo_ccqe_gr->SetPoint(j, ebo[j], 1);
	byEv_ca5_cc1pi_gr->SetPoint(j, ca5[j], (weights->At(ca5nom) > 0 ? weights->At(47+j)/weights->At(ca5nom) : 1));
	byEv_ca5_ncpiz_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpipm_gr->SetPoint(j, ca5[j], 1);
	byEv_manff_cc1pi_gr->SetPoint(j, manff[j], (weights->At(manffnom) > 0 ? weights->At(60+j)/weights->At(manffnom) : 1));
	byEv_manff_ncpiz_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpipm_gr->SetPoint(j, manff[j], 1);
	byEv_bgscl_cc1pi_gr->SetPoint(j, bgscl[j], (weights->At(bgsclnom) > 0 ? weights->At(73+j)/weights->At(bgsclnom) : 1));
        byEv_bgscl_ncpiz_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpipm_gr->SetPoint(j, bgscl[j], 1);
        byEv_dismpishp_ccoth_gr->SetPoint(j, dismpishp[j], 1);
      }
      for(int j=13; j<21; ++j) byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
    }

    else if(abs(mode) >= 31 && abs(mode) <= 32) { // ncpi0: ca5, manff, bgscl
      for(int j=0; j<13; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
	byEv_pfo_ccqe_gr->SetPoint(j, pfo[j], 1);
	byEv_ebo_ccqe_gr->SetPoint(j, ebo[j], 1);
	byEv_ca5_cc1pi_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpiz_gr->SetPoint(j, ca5[j], (weights->At(ca5nom) > 0 ? weights->At(47+j)/weights->At(ca5nom) : 1));
	byEv_ca5_ncpipm_gr->SetPoint(j, ca5[j], 1);
	byEv_manff_cc1pi_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpiz_gr->SetPoint(j, manff[j], (weights->At(manffnom) > 0 ? weights->At(60+j)/weights->At(manffnom) : 1));
	byEv_manff_ncpipm_gr->SetPoint(j, manff[j], 1);
	byEv_bgscl_cc1pi_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpiz_gr->SetPoint(j, bgscl[j], (weights->At(bgsclnom) > 0 ? weights->At(73+j)/weights->At(bgsclnom) : 1));
        byEv_bgscl_ncpipm_gr->SetPoint(j, bgscl[j], 1);
        byEv_dismpishp_ccoth_gr->SetPoint(j, dismpishp[j], 1);
      }
      for(int j=13; j<21; ++j) byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
    }

    else if(abs(mode) >= 33 && abs(mode) <= 34) { // ncpipm: ca5, manff, bgscl
      for(int j=0; j<13; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
	byEv_pfo_ccqe_gr->SetPoint(j, pfo[j], 1);
	byEv_ebo_ccqe_gr->SetPoint(j, ebo[j], 1);
	byEv_ca5_cc1pi_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpiz_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpipm_gr->SetPoint(j, ca5[j], (weights->At(ca5nom) > 0 ? weights->At(47+j)/weights->At(ca5nom) : 1));
	byEv_manff_cc1pi_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpiz_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpipm_gr->SetPoint(j, manff[j], (weights->At(manffnom) > 0 ? weights->At(60+j)/weights->At(manffnom) : 1));
	byEv_bgscl_cc1pi_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpiz_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpipm_gr->SetPoint(j, bgscl[j], (weights->At(bgsclnom) > 0 ? weights->At(73+j)/weights->At(bgsclnom) : 1));
        byEv_dismpishp_ccoth_gr->SetPoint(j, dismpishp[j], 1);
      }
      for(int j=13; j<21; ++j) byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
    }

    else if(abs(mode) >=17 && abs(mode) < 30) { // ccother: disshape
      for(int j=0; j<13; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
	byEv_pfo_ccqe_gr->SetPoint(j, pfo[j], 1);
	byEv_ebo_ccqe_gr->SetPoint(j, ebo[j], 1);
	byEv_ca5_cc1pi_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpiz_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpipm_gr->SetPoint(j, ca5[j], 1);
	byEv_manff_cc1pi_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpiz_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpipm_gr->SetPoint(j, manff[j], 1);
	byEv_bgscl_cc1pi_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpiz_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpipm_gr->SetPoint(j, bgscl[j], 1);
        byEv_dismpishp_ccoth_gr->SetPoint(j, dismpishp[j], (weights->At(dismpishpnom) > 0 ? weights->At(86+j)/weights->At(dismpishpnom) : 1));
      }
      for(int j=13; j<21; ++j) byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
    }

    else { // none for all other modes
      for(int j=0; j<13; ++j) {
	byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
	byEv_pfo_ccqe_gr->SetPoint(j, pfo[j], 1);
	byEv_ebo_ccqe_gr->SetPoint(j, ebo[j], 1);
	byEv_ca5_cc1pi_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpiz_gr->SetPoint(j, ca5[j], 1);
	byEv_ca5_ncpipm_gr->SetPoint(j, ca5[j], 1);
	byEv_manff_cc1pi_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpiz_gr->SetPoint(j, manff[j], 1);
	byEv_manff_ncpipm_gr->SetPoint(j, manff[j], 1);
	byEv_bgscl_cc1pi_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpiz_gr->SetPoint(j, bgscl[j], 1);
        byEv_bgscl_ncpipm_gr->SetPoint(j, bgscl[j], 1);
        byEv_dismpishp_ccoth_gr->SetPoint(j, dismpishp[j], 1);
      }
      for(int j=13; j<21; ++j) byEv_maqe_ccqe_gr->SetPoint(j, maqe[j], 1);
    }
    spltree->Fill();

  } // end of event loop

  outfile->cd();
  spltree->Write();
  outfile->Close();
  return 0;

}

int ParseArgs(int argc, char **argv)
{
  int nargs = 4;
  fUseWildcard = false;
  for(int i = 1; i <= 2*nargs; ++i) {
    if(std::string(argv[i]) == "-i") { fInputFileName = argv[i+1]; ++i; }
    else if(std::string(argv[i]) == "-s") {fInputSKFileName = argv[i+1]; ++i; }
    else if(std::string(argv[i]) == "-o") {fOutputFileName = argv[i+1]; ++i; }
    else if(std::string(argv[i]) == "-w") {fUseWildcard = atoi(argv[i+1]); ++i;}
    else {
      std::cout<<"Invalid argument: "<<argv[i]<<" "<<argv[i+1]<<std::endl;
      Usage();
    }
  }
  return 0;
}

void Usage()
{
  std::cout << "Cmd line syntax should be: \n"
	    << "genXsecSpline.exe -i [T2KReWeight-style wgt input file] -o outputfile -s [sk file] -w [0: process individual files; 1: produce a single file for all sk mc]"
	    << std::endl;
}
