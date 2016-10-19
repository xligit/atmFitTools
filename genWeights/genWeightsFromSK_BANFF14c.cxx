//
// 2015 NIWG-oscillation analysis dials for spline generation.
// See TN192 for details
//
//  Run with: 
//
// ./generateWeightsFromSK_NIWG2015.exe -s sk_inputfile [optional nevents] 
//
//
// Uses NIWGReWeight2015 and NEUT 3.3.3. Source the correct versions of 
// everything by doing:
// 
// Systematics for splines (in order):
//
// These ones have 13 dials (-3 sigma to 3 sigma in 0.5 sigma steps):
// MaQE
// Fermi momentum, pF (O)
// Binding energy, Eb (O)
// CA5
// 1pi axial form factor, MaNFFRES
// Isospin=1/2 background normalisation (BgSclRES)
// CC other shape uncertainty, dismpishp
// Second class currents (vector), SCCVecQE
// Second class currents (axial), SCCAxlQE
//
// This one has 1 dial only (tuning on):
// Rel RPA (note: calculate weights relative to untuned MC. Eventually will be applied on top of SF->RFG tuning and Non-rel RPA tuning, but we do not consider the errors associated with these two)
// 
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMath.h"

#include "T2KReWeight.h"
#include "T2KSyst.h"

#include "T2KGenieReWeight.h" 
#include "T2KNeutReWeight.h"
#include "T2KJNuBeamReWeight.h"
#include "T2KNIWGReWeight.h"

#include "SK__h1.h"
// For weight storer class
#include "T2KWeightsStorer.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2krew;

int fNskEvts = -1;
TString fSKFileName;
TString fInputParFileName;
TString fOutputFileName;

const int nPars = 7;
TString parNames[nPars] = {"MAQE","pF_O","EB_O","CA5","MANFFRES","BgRES","DISMPISHP"};
TString dialNames[nPars] = {"NXSec_MaCCQE", "NIWG2014a_pF_O16", "NIWG2014a_Eb_O16", "NXSec_CA5RES","NXSec_MaNFFRES","NXSec_BgSclRES","NIWG2012a_dismpishp"};
double sigma[nPars];
double nominal[nPars];
double prefit[nPars];
// set lower and upper bound by hand :(
double lb[nPars] = {0, 0.89, 0.44, 0, 0, 0, -9999};
double up[nPars] = {9999, 1.22, 1.56, 9999, 9999, 9999, 9999};

void Usage();
int ParseArgs(int argc, char *argv[]);

int main(int argc, char *argv[])
{
 
  // process the arguments
  int args = ParseArgs(argc, argv);
  if(args != 0){
    std::cerr << "Usage " << std::endl;
    return 0;
  }
  //T2KWeightsStorer *storer = new T2KWeightsStorer("testweights_skexample_nue.root"); //forweightstorer
  T2KWeightsStorer *storer = new T2KWeightsStorer(fOutputFileName);
  
  cout << "Starting to reweight NIWG events from SK file: " << fSKFileName << endl;
    
  // Load the SK "h1" tree
  TFile * SK_infile = new TFile(fSKFileName, "OPEN");
  if(!SK_infile){
    cerr << "Cannot open SK file!" << endl;
   exit(1);
  }
  TTree * SK_tree;
  SK_tree = (TTree*) SK_infile->Get("h1");
  if(!SK_tree){
    cerr << "Cannot find tree h1! Looking for mtuple instead" << endl;
    SK_tree = (TTree*) SK_infile->Get("mtuple");
  }
  if(!SK_tree){
    cerr << "Cannot find tree mtuple! Looking for minituple instead" << endl;
    SK_tree = (TTree*) SK_infile->Get("minituple");
  }
  if(!SK_tree){
    cerr << "Cannot find tree mtuple! This is a real problem." << endl;
  }  
  std::cout << "Number of entries in SK_tree: " << SK_tree->GetEntries() << std::endl;
  
  // Instantiate the reader of the SK tree class
  SK::SK__h1 *skVtxs = new SK::SK__h1(SK_tree,1);
  std::cout << "Number of entries in skVtxs: " << skVtxs->GetEntries() << std::endl;
  
  const int ndraws = 13*nPars+8;//13*9+2+7;//+21*2+16; // -3 sigma ... +3 sigma in 0.5sigma steps (13 variations) for 9 parameters, plus 2 RPA tunings, plus 7 extra tweaks (-7 sigma to -3.5 sigma) for MaQE, plus SF-mode MaQE and pFsf (both 10 sigma in 0.5 sigma steps), plus MEC q3 (16 steps)

  // input parameters file
  TFile *parfile = new TFile(fInputParFileName, "OPEN");
  if(!parfile->IsOpen()) exit(1);  
  cout << "Loading input parameters file: " << fInputParFileName << endl << endl;  
  // Get included parameters from file
  TObjArray *param_list_arr = (TObjArray*)parfile->Get("param_list");
  TMatrixDSym *prefit_cov = (TMatrixDSym*)parfile->Get("prefit_cov");
  TVectorD *prefit_params = (TVectorD*)parfile->Get("prefit_params");

  for (int i=0; i<param_list_arr->GetEntries(); i++) {
    TObjString *param_list_obj = (TObjString*)param_list_arr->At(i);
    TString paramStr = param_list_obj->GetString().Data();// name of the parameter in BANFF output file
    if (paramStr.Contains("SK")) continue;
    for (int ipar=0; ipar<nPars; ipar++) {
      if (parNames[ipar] == paramStr) {
	if (paramStr.Contains("DISMPISHP")) nominal[ipar] = 0;
	else nominal[ipar] = 1;
	sigma[ipar] = TMath::Sqrt((*prefit_cov)(i,i));
	prefit[ipar] = (*prefit_params)(i);
      }
    }
  }

    
#ifdef __T2KRW_NIWG_ENABLED__
#ifdef __T2KRW_NEUT_ENABLED__
  // Make a t2kreweighting object and add a NIWG and NEUT weighting engine. 
  t2krew::T2KReWeight rw; 
  rw.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());
  rw.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
  
  // NIWG 2015
  // Uncertainties
  rw.Systematics().Include(t2krew::kNXSec_MaCCQE);
  rw.Systematics().Include(t2krew::kNIWG2014a_pF_O16);
  rw.Systematics().Include(t2krew::kNIWG2014a_Eb_O16); 
  rw.Systematics().Include(t2krew::kNXSec_CA5RES);
  rw.Systematics().Include(t2krew::kNXSec_MaNFFRES);
  rw.Systematics().Include(t2krew::kNXSec_BgSclRES);
  rw.Systematics().Include(t2krew::kNIWG2012a_dismpishp);
  //rw.Systematics().Include(t2krew::kNXSec_SCCVecQE);
  //rw.Systematics().Include(t2krew::kNXSec_SCCAxlQE);
  //rw.Systematics().Include(t2krew::kNIWGSpectralFunc_pFsfO16_smooth);
  //rw.Systematics().Include(t2krew::kNIWGMEC_q3Cut);
  
  // Tunings
  //rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_norm);
  //rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_shape);
  rw.Systematics().Include(t2krew::kNXSec_VecFFCCQE); // used to set MAQE to act according to for RFG MC (2) or SF MC (402). Should always be set to 2 to ensure that we get the correct splines when SF -> RFG tuning is applied later, except for in the case of SF mode MAQE.  
  
  // Absolute tweak dials set the fractional uncertainty, instead of 
  // in units of "sigma", defined in the code.
  // Useful so that you define the uncertainty within the code, as what is
  // hardcoded may not be the same as what is used for analysis. 
  
  // Uncertanties:
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_MaCCQE);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_pF_O16);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_Eb_O16); 
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_CA5RES);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_MaNFFRES);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_BgSclRES);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_dismpishp);
  //rw.Systematics().SetAbsTwk(t2krew::kNXSec_SCCVecQE);
  //rw.Systematics().SetAbsTwk(t2krew::kNXSec_SCCAxlQE);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWGSpectralFunc_pFsfO16_smooth);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_q3Cut);
   
  // Tunings:
  //rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_norm);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_shape);
  
  std::cout << std::endl
            << "-------------------------------------------------" << std::endl
            << "1-sigma errors are set as: " << std::endl
            << "MaQE: " << sigma[0] << std::endl
            << "pF (O): put in by hand between boundaries" << sigma[1] << std::endl
            << "Eb (O): put in by hand between boundaries" << sigma[2] << std::endl
            << "Ca5: " << sigma[3] << std::endl
            << "MaNFFRES: " << sigma[4] << std::endl
            << "BgSclRES: " << sigma[5] << std::endl
            << "Dismpishp: " << sigma[6] << std::endl
            << "-------------------------------------------------"
            << std::endl << std::endl;

  // loop over total variations, changing each tweak dial accordingly
  for(int dial=0; dial<=ndraws; dial++){
  
    rw.Systematics().SetTwkDial(t2krew::kNXSec_VecFFCCQE, 2); // Note that when we set MAQE, we need to set NXSec_VecFFCCQE to 2 for SF->RFG MC. Should be set to 2 all the time, except for SF MAQE variations (127 - 147) to ensure we get the correct behaviour when SF->RFG tuning is applied later.
    double maqedial = 0;
    double pfodial = 0;
    double ebodial = 0;
    double cadial = 0;
    double manresdial = 0;
    double bkgdial = 0;
    double dismpidial = 0;

    // Set uncertainties
    if(dial>=0&&dial<=20){
      maqedial = (dial - 14)*0.5*sigma[0];
    }else if(dial>=21&&dial<=33){
      if (dial < 27) pfodial = (dial-21-6)*0.5*0.11/3.0;
      else           pfodial = (dial-21-6)*0.5*0.22/3.0;
    }else if(dial>=34&&dial<=46){
      ebodial = (dial-34-6)*0.5*0.56/3.0;
    }else if(dial>=47&&dial<=59){
      cadial = (dial-47-6)*0.5*sigma[3]; 
    }else if(dial>=60&&dial<=72){
      manresdial = (dial-60-6)*0.5*sigma[4]; 
    }else if(dial>=73&&dial<=85){
      bkgdial = (dial-73-6)*0.5*sigma[5]; 
    }else if(dial>=86&&dial<=98){
      dismpidial = (dial-86-6)*0.5*sigma[6]; 
    }

    // Set uncertainty dials
    rw.Systematics().SetTwkDial(t2krew::kNXSec_MaCCQE, maqedial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_pF_O16, pfodial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_Eb_O16, ebodial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_CA5RES, cadial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_MaNFFRES, manresdial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_BgSclRES, bkgdial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2012a_dismpishp, dismpidial);

    rw.Reconfigure();
    storer->NewSystSet(rw.Systematics()); // save the current twk dial values to the weight storer class

   // Loop over SK entries and calculate weights.
   if(fNskEvts < 0) fNskEvts = skVtxs->GetEntries();
   
   for(int i = 0; i < fNskEvts; i++){

   skVtxs->GetEntry(i);
   Double_t weight =1.0;
   weight = rw.CalcWeight(skVtxs);
   storer->AddWeight(weight); // add weight for each
   } // event loop

 } // index of tweak dial changes


  cout << "Outputfile testweights_skexample_nue.root has weight storer tree"  << endl;


 storer->SaveToFile(); // save the weights to a file
 #endif // __T2KRW_NIWG_ENABLED__
 #endif // __T2KRW_NEUT_ENABLED__
 delete storer; 
 
 return 0;
}

// Print the cmd line syntax
void Usage(){
  cout << "Cmd line syntax should be:" << endl;
  cout << "generateWeightsFromSK_NIWGexample.exe -s sk_inputfile [-e nevents] -i banff_file -o outfile" << endl;
  exit(1);
}

int ParseArgs(int argc, char **argv){
  int nargs = 3;
  if (argc < (nargs*2+1)) { Usage(); }
  for(int i = 1; i < argc; i++){
    if(string(argv[i]) == "-i") {fInputParFileName = argv[i+1]; i++;}
    else if(string(argv[i]) == "-s") {fSKFileName = argv[i+1]; i++; }
    else if(string(argv[i]) == "-e") {fNskEvts = std::atoi(argv[i+1]); i++;}
    else if(string(argv[i]) == "-o") {fOutputFileName = argv[i+1]; i++; }
    else {
      cout << "Invalid argument:" << argv[i] << " "<< argv[i+1] << endl;
      Usage();
    }
  }
  return 0;
}

