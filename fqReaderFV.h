//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun  9 15:12:53 2015 by ROOT version 5.28/00c
// from TTree h1fv/h1fv
// found on file: ccqe_clean.root
//////////////////////////////////////////////////////////
#ifndef __FQREADERFV_H__
#define __FQREADERFV_H__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TTree.h"
#include "TBranch.h"

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
