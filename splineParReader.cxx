#ifndef splineParReader_cxx
#define splineParReader_cxx


#include "splineParReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
 
///////////////////////////////////////////////
// Draw 1D plot of the spline for a given bin
TH1D* splineParReader::drawSpline(int ibin){
  TH1D* h  = new TH1D();
  return h;
}


/////////////////////////////////////////////
// Draw a 2D plot of the histogram variation
TH2D* splineParReader::draw2D(int ievent){
 
  if (ievent>0) fChain->GetEntry(ievent);

  // setup histogram
  int nbinsx = nhistobins;
  double xmin   = 0;
  double xmax   = nhistobins;
  int nbinsy = npoints;
  double ymin = systParValues[0];
  double ymax = systParValues[npoints-1];
  TH2D* h2 = new TH2D("h2d","h2d",nbinsx,xmin,xmax,nbinsy,ymin,ymax);

  // fill histogram
  for (int iy=0; iy<npoints; iy++){
    for (int ix=1; ix<=nbinsx; ix++){
       double value = binWeight[iy][ix];
       h2->SetBinContent(ix,iy+1,value);
    }
  }

  // draw histogram
  h2->Draw("lego2");

  return h2;
}

void splineParReader::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L splineParReader.cxx
//      Root > splineParReader t
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

#endif
