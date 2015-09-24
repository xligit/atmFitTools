//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 24 15:54:13 2015 by ROOT version 5.28/00c
// from TTree MCMCpath/MCMCpath
// found on file: mcmc1000nom.root
//////////////////////////////////////////////////////////

#ifndef deletethis_h
#define deletethis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class deletethis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           npars;
   Int_t           step;
   Float_t         par[100];

   // List of branches
   TBranch        *b_npars;   //!
   TBranch        *b_step;   //!
   TBranch        *b_pars;   //!

   deletethis(TTree *tree=0);
   virtual ~deletethis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef deletethis_cxx
deletethis::deletethis(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcmc1000nom.root");
      if (!f) {
         f = new TFile("mcmc1000nom.root");
      }
      tree = (TTree*)gDirectory->Get("MCMCpath");

   }
   Init(tree);
}

deletethis::~deletethis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t deletethis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t deletethis::LoadTree(Long64_t entry)
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

void deletethis::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npars", &npars, &b_npars);
   fChain->SetBranchAddress("step", &step, &b_step);
   fChain->SetBranchAddress("par", par, &b_pars);
   Notify();
}

Bool_t deletethis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void deletethis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t deletethis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef deletethis_cxx
