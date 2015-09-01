//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun  9 15:12:53 2015 by ROOT version 5.28/00c
// from TTree h1fv/h1fv
// found on file: ccqe_clean.root
//////////////////////////////////////////////////////////

#ifndef fqReaderFV_h
#define fqReaderFV_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class fqReaderFV {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           mode;
   UInt_t          nhitac;
   Int_t           absmode;
   Int_t           npar;
   Int_t           wallv;
   Int_t           ipv[50];
   Float_t         posv[3];
   Float_t         dirv[50][3];
   Float_t         pmomv[50];
   Float_t         vwall[50];
   Float_t         vtowall[50];
   Float_t         vphiwall[50];
   Float_t         fqitwnd[10];
   Float_t         fqtotq[10];
   Float_t         fq1rpcflg[10][7];
   Float_t         fq1rmom[10][7];
   Float_t         fq1rt0[10][7];
   Float_t         fq1rtotmu[10][7];
   Float_t         fq1rnll[10][7];
   Float_t         fq1rpos[10][7][3];
   Float_t         fq1rdir[10][7][3];
   Float_t         fqmsnfit;
   Float_t         fqmspcflg[5];
   Float_t         fqmstotmu[5];
   Float_t         fqmsnll[5];
   Float_t         fqmsmom[5][20];
   Float_t         fqmspos[5][20];
   Float_t         fqmsdir[5][20];
   Int_t           fqnse;
   Float_t         fq1rwall[10][7];
   Float_t         fq1rtowall[10][7];
   Float_t         fq1rphiwall[10][7];
   Float_t         fq1rcornerangle[10][7];
   Int_t           fq1rentry[10][7];
   Float_t         fqmrnring[100];
   Int_t           vdecayenpar;
   Float_t         vdecayemom[10];
   Float_t         vdecayepos[10][3];
   Float_t         vdecayedir[10][3];
   Float_t         vdecayewall[10];
   Float_t         vdecayetowall[10];
   Float_t         vdecayephiwall[10];

   // List of branches
   TBranch        *b_mode;   //!
   TBranch        *b_nhitac;   //!
   TBranch        *b_absmode;   //!
   TBranch        *b_npar;   //!
   TBranch        *b_wallv;   //!
   TBranch        *b_ipv;   //!
   TBranch        *b_posv;   //!
   TBranch        *b_dirv;   //!
   TBranch        *b_pmomv;   //!
   TBranch        *b_vwall;   //!
   TBranch        *b_vtowall;   //!
   TBranch        *b_vphiwall;   //!
   TBranch        *b_fqitwnd;   //!
   TBranch        *b_fqtotq;   //!
   TBranch        *b_fq1rpcflg;   //!
   TBranch        *b_fq1rmom;   //!
   TBranch        *b_fq1rt0;   //!
   TBranch        *b_fq1rtotmu;   //!
   TBranch        *b_fq1rnll;   //!
   TBranch        *b_fq1rpos;   //!
   TBranch        *b_fq1rdir;   //!
   TBranch        *b_fqmsnfit;   //!
   TBranch        *b_fqmspcflg;   //!
   TBranch        *b_fqmstotmu;   //!
   TBranch        *b_fqmsnll;   //!
   TBranch        *b_fqmsmom;   //!
   TBranch        *b_fqmspos;   //!
   TBranch        *b_fqmsdir;   //!
   TBranch        *b_fqnse;   //!
   TBranch        *b_fq1rwall;   //!
   TBranch        *b_fq1rtowall;   //!
   TBranch        *b_fq1rphiwall;   //!
   TBranch        *b_fq1rcornerangle;   //!
   TBranch        *b_fq1rentry;   //!
   TBranch        *b_fqmrnring;   //!
   TBranch        *b_vdecayenpar;   //!
   TBranch        *b_vdecayemom;   //!
   TBranch        *b_vdecayepos;   //!
   TBranch        *b_vdecayedir;   //!
   TBranch        *b_vdecayewall;   //!
   TBranch        *b_vdecayetowall;   //!
   TBranch        *b_vdecayephiwall;   //!

   fqReaderFV(TTree *tree=0);
   virtual ~fqReaderFV();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fqReaderFV_cxx
fqReaderFV::fqReaderFV(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ccqe_clean.root");
      if (!f) {
         f = new TFile("ccqe_clean.root");
      }
      tree = (TTree*)gDirectory->Get("h1fv");

   }
   Init(tree);
}

fqReaderFV::~fqReaderFV()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fqReaderFV::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fqReaderFV::LoadTree(Long64_t entry)
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

void fqReaderFV::Init(TTree *tree)
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

   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("nhitac", &nhitac, &b_nhitac);
   fChain->SetBranchAddress("absmode", &absmode, &b_absmode);
   fChain->SetBranchAddress("npar", &npar, &b_npar);
   fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
   fChain->SetBranchAddress("ipv", ipv, &b_ipv);
   fChain->SetBranchAddress("posv", posv, &b_posv);
   fChain->SetBranchAddress("dirv", dirv, &b_dirv);
   fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
   fChain->SetBranchAddress("vwall", vwall, &b_vwall);
   fChain->SetBranchAddress("vtowall", vtowall, &b_vtowall);
   fChain->SetBranchAddress("vphiwall", vphiwall, &b_vphiwall);
   fChain->SetBranchAddress("fqitwnd", fqitwnd, &b_fqitwnd);
   fChain->SetBranchAddress("fqtotq", fqtotq, &b_fqtotq);
   fChain->SetBranchAddress("fq1rpcflg", fq1rpcflg, &b_fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", fq1rmom, &b_fq1rmom);
   fChain->SetBranchAddress("fq1rt0", fq1rt0, &b_fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", fq1rtotmu, &b_fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll[10][7]", fq1rnll, &b_fq1rnll);
   fChain->SetBranchAddress("fq1rpos", fq1rpos, &b_fq1rpos);
   fChain->SetBranchAddress("fq1rdir", fq1rdir, &b_fq1rdir);
   fChain->SetBranchAddress("fqmsnfit", &fqmsnfit, &b_fqmsnfit);
   fChain->SetBranchAddress("fqmspcflg", fqmspcflg, &b_fqmspcflg);
   fChain->SetBranchAddress("fqmstotmu", fqmstotmu, &b_fqmstotmu);
   fChain->SetBranchAddress("fqmsnll", fqmsnll, &b_fqmsnll);
   fChain->SetBranchAddress("fqmsmom", fqmsmom, &b_fqmsmom);
   fChain->SetBranchAddress("fqmspos", fqmspos, &b_fqmspos);
   fChain->SetBranchAddress("fqmsdir", fqmsdir, &b_fqmsdir);
   fChain->SetBranchAddress("fqnse", &fqnse, &b_fqnse);
   fChain->SetBranchAddress("fq1rwall", fq1rwall, &b_fq1rwall);
   fChain->SetBranchAddress("fq1rtowall", fq1rtowall, &b_fq1rtowall);
   fChain->SetBranchAddress("fq1rphiwall", fq1rphiwall, &b_fq1rphiwall);
   fChain->SetBranchAddress("fq1rcornerangle", fq1rcornerangle, &b_fq1rcornerangle);
   fChain->SetBranchAddress("fq1rentry", fq1rentry, &b_fq1rentry);
   fChain->SetBranchAddress("fqmrnring", fqmrnring, &b_fqmrnring);
   fChain->SetBranchAddress("vdecayenpar", &vdecayenpar, &b_vdecayenpar);
   fChain->SetBranchAddress("vdecayemom", vdecayemom, &b_vdecayemom);
   fChain->SetBranchAddress("vdecayepos", vdecayepos, &b_vdecayepos);
   fChain->SetBranchAddress("vdecayedir", vdecayedir, &b_vdecayedir);
   fChain->SetBranchAddress("vdecayewall", vdecayewall, &b_vdecayewall);
   fChain->SetBranchAddress("vdecayetowall", vdecayetowall, &b_vdecayetowall);
   fChain->SetBranchAddress("vdecayephiwall", vdecayephiwall, &b_vdecayephiwall);
   Notify();
}

Bool_t fqReaderFV::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fqReaderFV::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fqReaderFV::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fqReaderFV_cxx
