//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  8 18:55:13 2015 by ROOT version 5.28/00c
// from TTree splinePars/spinePars
// found on file: nom3_splineFactoryOut.root
//////////////////////////////////////////////////////////

#ifndef splineParReader_h
#define splineParReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "shared.h"

class splineParReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nbin;
   Int_t           nhistobins;
   Int_t           ncomponent;
   Int_t           nattribute;
   Int_t           nsample;
   Int_t           nmode;
   Int_t           nsystpar;
   Int_t           nsyspartot;
   Int_t           npoints;
   Double_t        systParValues[NPTSMAX];
   Double_t        binWeight[NPTSMAX][NHBINSMAX];

   // List of branches
   TBranch        *b_nbin;   //!
   TBranch        *b_nhistobins;   //!
   TBranch        *b_ncomponent;   //!
   TBranch        *b_nattribute;   //!
   TBranch        *b_nsample;   //!
   TBranch        *b_nmode;
   TBranch        *b_nsystpar;   //!
   TBranch        *b_nsyspartot;   //!
   TBranch        *b_npoints;   //!
   TBranch        *b_systParValues;   //!
   TBranch        *b_binWeight;   //!

   splineParReader(TTree *tree=0);
   virtual ~splineParReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef splineParReader_cxx
splineParReader::splineParReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("nom3_splineFactoryOut.root");
      if (!f) {
         f = new TFile("nom3_splineFactoryOut.root");
      }
      tree = (TTree*)gDirectory->Get("splinePars");

   }
   Init(tree);
}

splineParReader::~splineParReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t splineParReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t splineParReader::LoadTree(Long64_t entry)
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

void splineParReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("nbin", &nbin, &b_nbin);
   fChain->SetBranchAddress("nhistobins", &nhistobins, &b_nhistobins);
   fChain->SetBranchAddress("ncomponent", &ncomponent, &b_ncomponent);
   fChain->SetBranchAddress("nattribute", &nattribute, &b_nattribute);
   fChain->SetBranchAddress("nsample", &nsample, &b_nsample);
   if (fChain->GetListOfBranches()->FindObject("nmode"))
     fChain->SetBranchAddress("nmode", &nmode, &b_nmode);
   fChain->SetBranchAddress("nsystpar", &nsystpar, &b_nsystpar);
   fChain->SetBranchAddress("nsyspartot", &nsyspartot, &b_nsyspartot);
   fChain->SetBranchAddress("npoints", &npoints, &b_npoints);
   fChain->SetBranchAddress("systParValues", systParValues, &b_systParValues);
   fChain->SetBranchAddress("binWeight", binWeight, &b_binWeight);
   Notify();
}

Bool_t splineParReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void splineParReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t splineParReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef splineParReader_cxx
