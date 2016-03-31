//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 22 02:52:03 2016 by ROOT version 5.28/00h
// from TTree SplinesByEvent/SplinesByEvent
// found on file: /disk2/usr5/xiaoyue/skmc/test_fc/splines/jan14sk4_skdetsim13p90_neut532.reduc.000_fQv4r.splines.root
//////////////////////////////////////////////////////////

#ifndef splineReader_h
#define splineReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class splineReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
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
   TGraph          *byEv_sccvec_ccqe_gr;
   TGraph          *byEv_sccvec_ncoth_gr;
   TGraph          *byEv_sccaxl_ccqe_gr;
   TGraph          *byEv_sccaxl_ncoth_gr;
   TGraph          *byEv_rpa_ccqe_gr;
   
   splineReader(TTree *tree=0);
   virtual ~splineReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef splineReader_cxx
splineReader::splineReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/disk2/usr5/xiaoyue/skmc/test_fc/splines/jan14sk4_skdetsim13p90_neut532.reduc.000_fQv4r.splines.root");
      if (!f) {
         f = new TFile("/disk2/usr5/xiaoyue/skmc/test_fc/splines/jan14sk4_skdetsim13p90_neut532.reduc.000_fQv4r.splines.root");
      }
      tree = (TTree*)gDirectory->Get("SplinesByEvent");

   }
   Init(tree);
}

splineReader::~splineReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t splineReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t splineReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void splineReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetBranchAddress("byEv_maqe_ccqe_gr", &byEv_maqe_ccqe_gr);
   fChain->SetBranchAddress("byEv_pfo_ccqe_gr", &byEv_pfo_ccqe_gr); 
   fChain->SetBranchAddress("byEv_ebo_ccqe_gr", &byEv_ebo_ccqe_gr); 
   fChain->SetBranchAddress("byEv_ca5_cc1pi_gr", &byEv_ca5_cc1pi_gr); 
   fChain->SetBranchAddress("byEv_ca5_ncpiz_gr", &byEv_ca5_ncpiz_gr); 
   fChain->SetBranchAddress("byEv_ca5_ncpipm_gr", &byEv_ca5_ncpipm_gr); 
   fChain->SetBranchAddress("byEv_manff_cc1pi_gr", &byEv_manff_cc1pi_gr); 
   fChain->SetBranchAddress("byEv_manff_ncpiz_gr", &byEv_manff_ncpiz_gr); 
   fChain->SetBranchAddress("byEv_manff_ncpipm_gr", &byEv_manff_ncpipm_gr); 
   fChain->SetBranchAddress("byEv_bgscl_cc1pi_gr", &byEv_bgscl_cc1pi_gr); 
   fChain->SetBranchAddress("byEv_bgscl_ncpiz_gr", &byEv_bgscl_ncpiz_gr); 
   fChain->SetBranchAddress("byEv_bgscl_ncpipm_gr", &byEv_bgscl_ncpipm_gr); 
   fChain->SetBranchAddress("byEv_dismpishp_ccoth_gr", &byEv_dismpishp_ccoth_gr); 
   fChain->SetBranchAddress("byEv_sccvec_ccqe_gr", &byEv_sccvec_ccqe_gr); 
   fChain->SetBranchAddress("byEv_sccvec_ncoth_gr", &byEv_sccvec_ncoth_gr); 
   fChain->SetBranchAddress("byEv_sccaxl_ccqe_gr", &byEv_sccaxl_ccqe_gr); 
   fChain->SetBranchAddress("byEv_sccaxl_ncoth_gr", &byEv_sccaxl_ncoth_gr); 
   fChain->SetBranchAddress("byEv_rpa_ccqe_gr", &byEv_rpa_ccqe_gr); 

   Notify();
}

Bool_t splineReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void splineReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t splineReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef splineReader_cxx
