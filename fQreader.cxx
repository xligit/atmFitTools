#define fQreader_cxx
#include "fQreader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void fQreader::FillMap()
{
  std::cout<<"filling map"<<std::endl;
  fChain->GetEntry(0);
  double x,y;
  for (int i = 0; i < byEv_maqe_ccqe_gr->GetN(); ++i) {
    byEv_maqe_ccqe_gr->GetPoint(i,x,y);
    maqe[(double)i*0.5-7] = x;
    std::cout<<((double)i*0.5-7)<<"="<<maqe[(double)i*0.5-7]<<" ";
  }
  std::cout<<"maqe"<<std::endl;
  for (int i = 0; i < byEv_pfo_ccqe_gr->GetN(); ++i) {
    byEv_pfo_ccqe_gr->GetPoint(i,x,y);
    pf_o[(double)i*0.5-3.0] = x;
    std::cout<<((double)i*0.5-3.0)<<"="<<pf_o[(double)i*0.5-3.0]<<" ";
  }
  std::cout<<"pfo"<<std::endl;
  for (int i = 0; i < byEv_ebo_ccqe_gr->GetN(); ++i) {
    byEv_ebo_ccqe_gr->GetPoint(i,x,y);
    eb_o[(double)i*0.5-3.0] = x;
    std::cout<<((double)i*0.5-3.0)<<"="<<eb_o[(double)i*0.5-3.0]<<" ";
  }
  std::cout<<"ebo"<<std::endl;
  for (int i = 0; i < byEv_ca5_cc1pi_gr->GetN(); ++i) {
    byEv_ca5_cc1pi_gr->GetPoint(i,x,y);
    ca5[(double)i*0.5-3.0] = x;
    std::cout<<((double)i*0.5-3.0)<<"="<<ca5[(double)i*0.5-3.0]<<" ";
  }
  std::cout<<"ca5"<<std::endl;
  for (int i = 0; i < byEv_manff_cc1pi_gr->GetN(); ++i) {
    byEv_manff_cc1pi_gr->GetPoint(i,x,y);
    manffres[(double)i*0.5-3.0] = x;
    std::cout<<((double)i*0.5-3.0)<<"="<<manffres[(double)i*0.5-3.0]<<" ";
  }
  std::cout<<"manff"<<std::endl;
  for (int i = 0; i < byEv_bgscl_cc1pi_gr->GetN(); ++i) {
    byEv_bgscl_cc1pi_gr->GetPoint(i,x,y);
    bgres[(double)i*0.5-3.0] = x;
    std::cout<<((double)i*0.5-3.0)<<"="<<bgres[(double)i*0.5-3.0]<<" ";
  }
  std::cout<<"bgscl"<<std::endl;
  for (int i = 0; i < byEv_dismpishp_ccoth_gr->GetN(); ++i) {
    byEv_dismpishp_ccoth_gr->GetPoint(i,x,y);
    dismpishp[(double)i*0.5-3.0] = x;
    std::cout<<((double)i*0.5-3.0)<<"="<<dismpishp[(double)i*0.5-3.0]<<" ";
  }
  std::cout<<"dismpishp"<<std::endl;
  for (int i = 0; i < byEv_rpa_ccqe_gr->GetN(); ++i) {
    byEv_rpa_ccqe_gr->GetPoint(i,x,y);
    rpa[(double)i] = x;
  }
  std::cout<<"rpa"<<std::endl;
}

void fQreader::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L fQreader.C
//      Root > fQreader t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
