///////////////////////////////////////////////////////////////////////////////////
//
// Title: genWeights_2015.cxx
//
// Purpose: Generate weights for ND280 or SK sample, using the BANFF
//          parameters fitted central values (CVs) and covariance
//
// Usage:
//   ./genWeights_2015.exe -i <inputfile> -p <banff_parameter_file> -o <weight_outputfile>
//                          -horn <+1:nu-mode, -1:antinu-mode> -app <1 for appearance sample, 0 otherwise>
//                          -t <# of throws> -r <random seed for throws>
//                          -dslist <disable_sys_list>
//                          [--use-prefit --drop-flux --drop-xsec]
//
//      where,
//
//             <inputfile>: ND280 or SK MC file. Currently supported:
//                          1) NuMu CC (irods://QMULZone/home/asg/banff/numu/CC_MC_*.root)
//                          2) P0D NuE (http://nngroup.physics.sunysb.edu/~itaylor/T2K_Files/NuEAnalysisFiles/complete-Selection111024-[1-3].root
//                          3) Raw oaAnalysis file (only 1st vertex weight is stored)
//                          4) Raw SK tree (h1)
//                          5) NEUT neutroot file for e.g. MiniBooNE simulations (currently disabled because of issue with gcc4.6.3)
//
//  <banff_parameter_file>: Parameter file output by the BANFF fit
//
//           <# of throws>: Number of weights to generate for every event, each based
//                          on a throw from the fitted covariance matrix.
//                          Default is 0, only CV weights stored.
//
//     <weight_outputfile>: ROOT file containing TTree of weights for each event.
//                          First set in the weight array ( weights->At(0) ) is the
//                          fitted central value. The following sets in the array
//                          ( weights->At(i, i>0) ) correspond to each throw from
//                          the covariance matrix, if any (set by '-t'). Use these
//                          weights and evaluate the RMS in each bin of your distributions
//                          to build an error envelope.
//
//  <+1:nu-mode, -1:antinu-mode>: Set whether the MC file is generated based on the nu-mode or antinu-mode flux
//
//  <1 for appearance sample, 0 otherwise>: When processing SK nue/nueb appearance sample, set this flag to 1
//
//      <disable_sys_list>: Text file that lists the systematic parameter dials to disable. Nominal value will be used when calculating weights.
//
//      --use-prefit: Use parameters and covariance that are inputs to BANFF fit
//      --drop-flux: Don't apply any flux weights, even if flux parameters are defined
//      --drop-xsec: Don't apply any xsec weights, even if xsec parameters are defined
//
// History:
//
// 2015-04-13: Created based on genWeights_2013a - Shimpei Tobayama
//               Has only been tested on SK ntuples so far
//
///////////////////////////////////////////////////////////////////////////////////

//#define DEBUG

#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TAxis.h"

#include "T2KReWeight.h"
#include "T2KSyst.h"

#include "T2KGenieReWeight.h"
#include "T2KGenieUtils.h"

#include "T2KNeutReWeight.h"
#include "T2KNeutUtils.h"

#include "T2KNIWGReWeight.h"
#include "T2KNIWGUtils.h"

#include "T2KJNuBeamReWeight.h"

#ifdef __T2KRW_OAANALYSIS_ENABLED__
#include "ND__NRooTrackerVtx.h"
#endif

#include "T2KWeightsStorer.h"

#include "SK__h1.h"

#ifdef __T2KRW_NEUT_ENABLED__
//#define NEUTROOT
#ifdef NEUTROOT
#include "neutrootTreeSingleton.h"
#endif
#endif

using namespace std;

///////////////////////////////////////////////////////////////////////////////////

// Define Supported Analyses
const int nAnas = 5;
string anaTreeNames[nAnas] = {
  // Keep all RooTrackerVtx based files below here
  "SummaryTree/SummaryTree", // numuCC
  "AnalysisP0DNuE",          // P0DNuE
  "TruthDir/NRooTrackerVtx", // raw_oaAnalysis (this should be last in this section)
  // Keep all RooTrackerVtx based files above here
  
  // Keep all SK__h1 based files below here
  "h1",                      // raw_sk
  // Keep all SK__h1 based files above here
  
  // Keep all neutroot based files below here
  "neuttree"                 // raw_neutroot
  // Keep all neutroot based files above here
};

// This must correspond to the above "anaTreeNames"
enum anaEnum {numuCC,
  P0DNuE,
  raw_oaAnalysis,
  raw_sk,
  raw_neutroot
};
// Note: "raw_oaAnalysis" is last in the ND280 trees since some analyses have
//       the same "TruthDir/NRooTrackerVtx" tree structure

int iana=0; // Counter and flag for determining analysis

// Tree types T2KReWeight input
enum tree_type_enum {NRooTrackerVtx, // iana=0,1,2
  SK__h1,         // iana=3
  neutroot        // iana=4
};
int tree_type = -1;

// Parameters as Defined in the BANFF fit
const int nPars = 23;

// Name of the tuning parameters used in the BANFF output file.
// Prepend "Flx_" to flux parameters; detector type, horn mode and bin index will be appended in the code later
TString parNames[nPars] = {"Flx_Numu","Flx_Numub","Flx_Nue","Flx_Nueb","MAQE","pF_C","MEC_C","EB_C","pF_O","MEC_O","EB_O","MEC_NUBAR","CA5","MANFFRES","BgRES","SigmaNue","SigmaNuebar","DISMPISHP","CCCOH_C_0","CCCOH_O_0","NCCOH_0","NC1GAMMA_0","NCOTHER_FAR_0"};

// Name of the corresponding T2KReWeight dial for each of the parameters above
// For parameters whose value default to 0, not 1, make sure to add it in the "if" statement where the tweak dial is set later in the code
TString dialNames[nPars] = {"numu","numub","nue","nueb","NXSec_MaCCQE","NIWG2014a_pF_C12","NIWGMEC_Norm_C12","NIWG2014a_Eb_C12","NIWG2014a_pF_O16","NIWGMEC_Norm_O16","NIWG2014a_Eb_O16","NIWGMEC_Norm_other","NXSec_CA5RES","NXSec_MaNFFRES","NXSec_BgSclRES","NIWG2012a_ccnueE0","NIWG2012a_ccnueE0","NIWG2012a_dismpishp","NIWG2012a_cccohE0","NIWG2012a_cccohE0","NIWG2012a_nccohE0","NIWG2015b_ncgamma","NIWG2012a_ncotherE0"};

int iHornMode = 0;// +1:nu-mode, -1:antinu-mode

int iApp = 0;// 1: appearance sample, 0: otherwise

enum DetTypeEnum {Det_ND, Det_SK};
int iDetType = Det_SK;
TString jnuDet[2] = {"nd5","sk"};
TString strDetType[2]={"ND","SK"};

const int nNus = 4;
enum jnuFluxEnum {jnu_numu, jnu_numub, jnu_nue, jnu_nueb};
string jnuFlux[nNus] = {"numu", "numub", "nue", "nueb"};


// For determining if a parameter was included in BANFF fit and
// bin number in parameter vector/covariance matrix
int parIncluded[nPars] = {0};

// For checking binning against that hardcoded in reweighting engines
TAxis *parBins[nPars];

// Number of events in file to loop over (default is all events)
// Can be modified with '-e' command line argument
int fNEvts = -1;

// Filename containers
TString fInputFileName = "";
TString fOutputFileName = "";
TString fInputParsFileName = "";

TString fDisablSystList = "";

// Configuration of inputs
bool fUseFlux = true;
bool fUseXsec = true;
bool fUsePrefit = false;

// For matching to RooTrackerVtx
int anaNeutmode, NeutGenieIdx;
double anaEnu;

// For throwing parameters from covariance
int nThrows=0;
int rndSeed = 1867;

void Usage();
void ParseArgs(int argc, char **argv);

int main(int argc, char *argv[])
{
#if defined (__T2KRW_NEUT_ENABLED__) && defined (__T2KRW_JNUBEAM_ENABLED__) && defined (__T2KRW_NIWG_ENABLED__)
  ParseArgs(argc, argv);
  
  // Check that filenames were specified
  if (!fInputFileName.Length() || !fOutputFileName.Length() || !fInputParsFileName.Length()) Usage();
  
  // Check the MC mode flags
  if (!(iHornMode==1 || iHornMode==-1)) Usage();
  
  // Open the input file
  TFile * infile = new TFile(fInputFileName, "OPEN");
  if(!infile->IsOpen()) exit(1);
  
  cout << endl << "Loading input tree file: " << fInputFileName.Data() << endl << endl;
  
  TTree *input_tree = 0;
  for (iana=0; iana<nAnas; iana++) {
    input_tree = (TTree*) infile->Get(anaTreeNames[iana].c_str());
    if (input_tree) {
      cout << "Found " << anaTreeNames[iana] << " tree" << endl;
      break;
    }
  }
  
  if (!input_tree) {
    cerr << "Error: Cannot find input tree" << endl;
    exit (1);
  }
  
  // Add tree distinguishing here if implementing a new tree format
  if (iana<=raw_oaAnalysis) {
    tree_type = NRooTrackerVtx;
    iDetType = Det_ND;
    
#ifndef __T2KRW_OAANALYSIS_ENABLED__
    cout << "Error: oaAnalysis must be enabled upon build if using ND280 files" << endl;
    exit (1);
#endif
    
    cout << "Loading NRooTrackerVtx tree ";
  }
  else if (iana<=raw_sk) {
    tree_type = SK__h1;
    iDetType = Det_SK;
    cout << "Loading SK tree ";
  }
//  else if (iana<=raw_neutroot) {
//    tree_type = neutroot;
//    cout << "Loading neutroot tree ";
//  }
  else {
    cerr << "Error: Unknown analysis type, iana=" << iana << endl;
    exit (1);
  }
  
  // Specific analysis tree
  TTree *ana_tree = input_tree;
  
  // For matching to RooTrackerVtx
  TBranch *b_anaNeutmode, *b_anaEnu, *b_NeutGenieIdx;
  
  // Initialize numuCC analysis tree
  if (iana == numuCC) {
    input_tree = (TTree*)infile->Get("TruthDir/NRooTrackerVtx");
    
    ana_tree->SetBranchAddress("selNeutGenie", &NeutGenieIdx, &b_NeutGenieIdx);
    ana_tree->SetBranchAddress("ReactionCode", &anaNeutmode, &b_anaNeutmode);
    
  }
  // Initialize P0DNuE analysis tree
  else if (iana == P0DNuE) {
    input_tree = (TTree*)infile->Get("NRooTrackerVtx");
    
    ana_tree->SetBranchAddress("NEUTCode", &anaNeutmode, &b_anaNeutmode);
    ana_tree->SetBranchAddress("NuEnergy", &anaEnu,      &b_anaEnu);
  }
  
  if (!input_tree) {
    cout << "Error: Cannot find NRooTrackerVtx in analysis file" << endl;
    exit (1);
  }
  
  ana_tree->AddFriend(input_tree);
  
  // Initialize input tree (need to implement new tree formats here)
  int NVtx;
  TClonesArray * nRooVtxs;
  SK::SK__h1 *skVtxs;
  
  if (tree_type == NRooTrackerVtx) {
    nRooVtxs = new TClonesArray("ND::NRooTrackerVtx", 100);
    input_tree->SetBranchAddress("Vtx", &nRooVtxs);
    input_tree->SetBranchAddress("NVtx", &NVtx);
  }
  else if (tree_type == SK__h1) {
    skVtxs = new SK::SK__h1(input_tree,1);
    
    // Check for SK signal nue/nuebar file
    if (iApp) {
      cout << endl << "Setting flag as appearance sample..." << endl;
      skVtxs->filetype = 1;
    }
  }
#ifdef NEUTROOT
  else if (tree_type == neutroot) {
    // Create global tree object
    NeutrootTreeSingleton * neutrootTreeObj = NeutrootTreeSingleton::Instance(fInputFileName.Data());
  }
#endif
  
  // Get number of events in tree
  if(fNEvts < 0) fNEvts = input_tree->GetEntries();
  if (tree_type!=neutroot)
    cout << " with " << fNEvts << " events." << endl << endl;
  
  if (input_tree->GetEntries() != ana_tree->GetEntries()) {
    cout << "Error: Number of events in input_tree != ana_tree (" << ana_tree->GetEntries() << ")" << endl;
    exit (1);
  }
  
#ifdef __T2KRW_OAANALYSIS_ENABLED__
  vector<ND::NRooTrackerVtx*> allRooTrackerVtxs;
#endif
  
  vector<SK::SK__h1*> allSKVtxs;
  
  //====================================================================
  // Loop over events to store in memory
  for(int i = 0; i < fNEvts; i++){
    
    // Currently do not store neutroot events since it's difficult
    // in the current tree implementation
    if (tree_type == neutroot) break;
    
    ana_tree->GetEntry(i);
    
    // RooTrackerVtx
    if (tree_type == NRooTrackerVtx) {
#ifdef __T2KRW_OAANALYSIS_ENABLED__
      // Count number of matched vertices
      // (only useful if matching index "NeutGenieIdx" is not used)
      int nMatchedVertices=0;
      
      ND::NRooTrackerVtx * vtx=0;
      // Check if matching index was set properly
      if (NeutGenieIdx<0) {
        cout << "Warning: No matched RooTrackerVtx for event " << i << ", NeutGenieIdx = " << NeutGenieIdx << endl;
        
        // Temporarily assume there was a matched vertex but assign it a weight of 1
        nMatchedVertices=1;
        
      } else {
        
        // Loop over RooTrackerVtx's
        for(int j = 0; j<NVtx; j++){
          
          ND::NRooTrackerVtx* thisVtx = (ND::NRooTrackerVtx*) nRooVtxs->At(j);
          if(!thisVtx){
            cout << "Error: Cannot find NRooTrackerVtx object for event " << i << ", vertex " << j << endl;
            exit (1);
          }
          
          // For P0DNuE analysis, match vertex using NEUT mode and Enu
          if (iana==P0DNuE &&
              !(atoi(thisVtx->EvtCode->GetString().Data()) == anaNeutmode &&
                thisVtx->StdHepP4[0][3] == anaEnu))
            continue;
          
          // For numu CC analysis, use index variable
          else if ( iana==numuCC && j!=NeutGenieIdx)
            continue;
          
          vtx=thisVtx;
          
          nMatchedVertices++;
          
          // Only weight first vertex in raw oaAnalysis file
          // Since T2KWeightStorer tree is flat
          if (iana==raw_oaAnalysis) break;
        }
      }
      
      // Check number of matched vertices
      if (nMatchedVertices < 1) {
        cout << "Error: No matched RooTrackerVtx found in event " << i << endl;
        exit (1);
      }
      else if (nMatchedVertices > 1) {
        cout << "Error: " << nMatchedVertices << " matched RooTrackerVtx found in event " << i << endl;
        exit (1);
      }
      
      // Need to clone the vertex because TTree owns vtx (due to SetBranchAddress)
      if (vtx)
        allRooTrackerVtxs.push_back((ND::NRooTrackerVtx*)vtx->Clone());
      else
        allRooTrackerVtxs.push_back(0);
#endif
      
    } // end if tree_type==NRooTrackerVtx
    
    // SK tree
    else if (tree_type == SK__h1) {
      allSKVtxs.push_back(new SK::SK__h1(*skVtxs));
    }
  } // end loop over events
  
  //================================================================================
  
  // Initialize output weights tree and file
  t2krew::T2KWeightsStorer *wt_storer = new t2krew::T2KWeightsStorer(fOutputFileName.Data());
  
  if (iana==raw_oaAnalysis)
    cout << "Warning: Only first vertex in each event of raw oaAnalysis will be weighted" << endl << endl;
  
  
  // Load input parameters file
  TFile *parfile = new TFile(fInputParsFileName, "OPEN");
  if(!parfile->IsOpen()) exit(1);
  
  cout << "Loading input parameters file: " << fInputParsFileName << endl << endl;
  
  // Get included parameters from file
  TObjArray *param_list_arr = (TObjArray*)parfile->Get("param_list");
  for (int i=0; i<param_list_arr->GetEntries(); i++) {
    TObjString *param_list_obj = (TObjString*)param_list_arr->At(i);
    TString paramStr = param_list_obj->GetString().Data();// name of the parameter in BANFF output file
    
//    cout << "Enabling parameter: " << paramStr << endl;
    int parFound = 0;
    for (int ipar=0; ipar<nPars; ipar++) {
      TString tmpStr = parNames[ipar];
      if (tmpStr.Index("Flx_")==0) {// flux parameter
        tmpStr.ReplaceAll("Flx_","");
        
        TString strHornMode[2]={"NuMode","ANuMode"};
        for (int iDet=0; iDet<2; iDet++) {
          for (int iHrn=0; iHrn<2; iHrn++) {
            TString tmpFlxName = Form("%s%s%s",strDetType[iDet].Data(),strHornMode[iHrn].Data(),tmpStr.Data());
            if (paramStr.Index(tmpFlxName)==0) {// match!
              if (iDet==iDetType && iHrn==(iHornMode==-1 ? 1 : 0)) {// Detector and horn mode matched!
                if (parIncluded[ipar]==0) {// found first bin for this flux type!
                  parIncluded[ipar] = i+1;
                  parFound=1;
                  
                  TString tmpBinName=Form("%s_%s_%s_bins",jnuDet[iDetType].Data(),strHornMode[iHrn].Data(),tmpStr.Data());
                  tmpBinName.ToLower();
                  parBins[ipar] = (TAxis*)parfile->Get(tmpBinName.Data());
                  
                  if (!parBins[ipar]) {
                    cout << "Error: " << parNames[ipar] << " does not have associated TAxis in file." << endl;
                    exit (1);
                  }
                  else {
                    cout << Form("Enabling flux: %s",tmpBinName.Data()) << endl;
                  }
                }
                else parFound=-1;
              }
              else {// ignore flux for wrong detector & horn mode!
                parFound=-1;
              }
            }
          }
        }
      }
      else {// non-flux parameter
        int idxBinFlg = parNames[ipar].Index("Bins");
        if (idxBinFlg>0) {// binned parameter
          tmpStr = parNames[ipar](0,idxBinFlg);
        }
        
        if (paramStr.Index(tmpStr)==0) {// matched!
          if (parIncluded[ipar]==0) {
            parIncluded[ipar] = i+1;
            parFound = 1;
            
            if (idxBinFlg>0) {
              
              parBins[ipar] = (TAxis*)parfile->Get(Form("%sbins",tmpStr.Data()));
              
              if (!parBins[ipar]) {
                cout << "Error: " << parNames[ipar] << " does not have associated TAxis in file." << endl;
                exit (1);
              }
            }
          }
          else parFound=-1;
        }
      }
      
    }
    
    if (parFound==0) {
      cout << "Error: " << paramStr.Data() << " has not been implemented in this executable" << endl;
      exit (1);
    }
  }
  cout << endl;
  
  // disable dials on the list
  ifstream fList(fDisablSystList.Data());
  if (!fList.fail()) {
    
    string strLine;
    while (getline(fList,strLine)) {
      for (int ipar=0; ipar<nPars; ipar++) {
        if (parNames[ipar].CompareTo(strLine.c_str())==0) {
          cout << Form("Disabling dial for %s",parNames[ipar].Data()) << endl;
          parIncluded[ipar] = 0;
          parBins[ipar] = 0;
        }
      }
    }
    fList.close();
    cout << endl;
  }
  
  // Initialize T2KReWeight and engines
  t2krew::T2KReWeight rw;
  rw.AdoptWghtEngine("jnubeam_rw", new t2krew::T2KJNuBeamReWeight());
  rw.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());
  rw.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
  
  jnubeam::rew::JReWeightEnu2013a* jrw_ptr = (jnubeam::rew::JReWeightEnu2013a*)((t2krew::T2KJNuBeamReWeight*)rw.WghtEngine("jnubeam_rw"))->GetWghtCalc("enu2013a");
  if (!jrw_ptr) {
    cout << endl << "Error: \"enu2013a\" WghtCalc not enabled in \"jnubeam_rw\" WghtEngine" << endl;
    exit (-1);
  }
  
  // WARNING: There is no check of neutrino type, detector nor mode ordering
  // (assumed order is as given by parNames[] defined above)
  for (int ipar=0; ipar<nPars; ipar++) {
    if (!parBins[ipar]) continue;
    
    TAxis *axis;
    if (parNames[ipar].Index("Flx_")==0) {// flux parameter
      
      // Get neutrino type
      int jnuFluxType=-1;
      for (int inu=0; inu<nNus; inu++)
        if (dialNames[ipar].CompareTo(jnuFlux[inu].c_str(),TString::kIgnoreCase)==0) jnuFluxType = inu;
      
      if (jnuFluxType<0) {
        cout << "Error: Can't find neutrino type in " << parNames[ipar].Data() << endl;
        exit (1);
      }
      
      int FlxIdx = ( iDetType==Det_ND ? jnuFluxType : jnuFluxType+4 );
      jrw_ptr->axis[FlxIdx] = axis = parBins[ipar]; // Copy axis to jrw instead, so checks below are redundant
    }
    else {
      cout << "Histogrammed niwg parameters are not supported yet. Please manually add the parameter for each bin to the systematic list." << endl;
      cout << "Ensure the binning in the input file is consistent with what's hardcoded in NIWGReWeight!" << endl;
      exit (1);
    }
  }
  
  // Get fitted parameter central values
  TVectorD postfit_params(*(TVectorD*)parfile->Get(Form("%s_params",(fUsePrefit ? "prefit":"postfit"))));
  
  // Get fitted parameter covariance matrix
  TMatrixDSym *postfit_cov = (TMatrixDSym*)(parfile->Get(Form("%s_cov",(fUsePrefit ? "prefit":"postfit"))))->Clone();
  
  // Throwing class
  ThrowParms *par_throws = 0;
  if (nThrows>0) {
    par_throws = new ThrowParms(postfit_params, (*postfit_cov));
    par_throws->SetSeed(rndSeed);
    
    cout << "Starting " << nThrows << " throws with random seed = " << rndSeed << endl << endl;
  }
  
  
  for (int ithrow=0; ithrow<nThrows+1; ithrow++) {
    
    TVectorD curr_parset = postfit_params;
    
    // First parameter set is the fitted central value
    // Else get a random throw based on covariance
    if (ithrow!=0) {
      std::vector<double> par_throw;
      par_throws->ThrowSet(par_throw);
      for(int j=0; j<(int)par_throw.size(); j++)
        curr_parset(j) = par_throw[j];
    }
    
    // Enable and set T2KSyst values
    int nParsXnBins=0;
    for (int ipar=0; ipar<nPars; ipar++) {
      
      if (!parIncluded[ipar]) continue;
      
      nParsXnBins = parIncluded[ipar]-1;// For histogrammed parameters, this is the parameter index of the first bin
      
      t2krew::T2KSyst_t t2ksyst = t2krew::kSystNull;
      
      // Flux parameters
      if (parNames[ipar].Index("Flx_")==0) {// flux parameter
        
        // Get neutrino type
        int jnuFluxType=-1;
        for (int inu=0; inu<nNus; inu++)
          if (dialNames[ipar].CompareTo(jnuFlux[inu].c_str(),TString::kIgnoreCase)==0) jnuFluxType = inu;
        
        if (jnuFluxType<0) {
          cout << "Error: Can't find neutrino type in " << parNames[ipar].Data() << endl;
          exit (1);
        }
        
        // Enable T2KSyst's
        for (int ibin=0; ibin<parBins[ipar]->GetNbins(); ibin++) {
          
          // See ${T2KREWEIGHT}/src/T2KSyst.h for naming convention
          string systString = Form("JEnu2013a_%s%s%d",jnuDet[iDetType].Data(),jnuFlux[jnuFluxType].c_str(),ibin);
          
          t2ksyst = t2krew::T2KSyst::FromString(systString);
          
          if (t2ksyst == t2krew::kSystNull) {
            cout << "Error: Unknown T2KSyst: " << systString.c_str() << endl;
            exit (1);
          }
          
          //Only if using flux parameters
          if(fUseFlux){
            rw.Systematics().Include(t2ksyst);
            rw.Systematics().SetAbsTwk(t2ksyst);
            rw.Systematics().SetTwkDial(t2ksyst,curr_parset[nParsXnBins] - 1);
          }
          
          nParsXnBins++;
          
        }
      } // End flux parameters
      
      // NEUT & NIWG parameters
      else {
        
        // Enu binned parameters
        if (parBins[ipar]) {
          cout << "Binned NIWG parameters not supported yet!" << endl;
          exit (1);
        }
        // Non-Enu binned parameters
        else {
          
          string systString = dialNames[ipar].Data();
          t2ksyst = t2krew::T2KSyst::FromString(systString);
          
          if (t2ksyst == t2krew::kSystNull) {
            cout << "Error: Unknown T2KSyst: " << systString.c_str() << endl;
            exit (1);
          }
          
          //Only if using xsec parameters
          if(fUseXsec){
            rw.Systematics().Include(t2ksyst);
            rw.Systematics().SetAbsTwk(t2ksyst);
            
            // Parameters passed straight to NIWGReWeight
            if (parNames[ipar].CompareTo("DISMPISHP")==0)
              rw.Systematics().SetTwkDial(t2ksyst,curr_parset[nParsXnBins]);// parameter which defaults to 0
            else
              rw.Systematics().SetTwkDial(t2ksyst,curr_parset[nParsXnBins] - 1);// parameter which defaults to 1
          }
          
        }
        
      } // End NIWG Parameters
      
    }
    
    rw.Systematics().Include(t2krew::kNIWG2014a_SF_RFG);
    rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_norm);
    rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_shape);
    rw.Systematics().Include(t2krew::kNXSec_VecFFCCQE); // used to set MAQE to act according to for RFG MC (2) or SF MC (402). Should always be set to 2 to ensure that we get the correct splines when SF -> RFG tuning is applied later, except for in the case of SF mode MAQE.
    
    // RPA tuning takes AbsTwk too
    rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_norm);
    rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_shape);
    
    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_SF_RFG,1);
    rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm,1);
    rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,0);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_VecFFCCQE, 2); // Note that when we set MAQE, we need to set NXSec_VecFFCCQE to 2 for SF->RFG MC. Should be set to 2 all the time, except for SF MAQE variations (127 - 147) to ensure we get the correct behaviour when SF->RFG tuning is applied later.
    
    
    rw.Reconfigure();
    
    // Store the values of the current set of parameters
    wt_storer->NewSystSet(rw.Systematics());
    
    
    //////////////////////////////////////////////////////////
    
    
    // Loop over events in memory
    for(int i = 0; i < fNEvts; i++){
      
      //if (i%10000==0) cout << "Event " << i << endl;
      
      Double_t weight = 1.0;
      
      if (tree_type == NRooTrackerVtx) {
#ifdef __T2KRW_OAANALYSIS_ENABLED__
        ND::NRooTrackerVtx* vtx=allRooTrackerVtxs[i];
        if(vtx) weight = rw.CalcWeight(vtx);
#endif
      }
      
      else if (tree_type == SK__h1) {
        SK::SK__h1* vtx=allSKVtxs[i];
        if(vtx) weight = rw.CalcWeight(vtx);
      }
      
      else if (tree_type == neutroot) {
        weight = rw.CalcWeight(input_tree, i);
      }
      
      // Store weight for the event
      wt_storer->AddWeight(weight);
      
    } // end loop over events
    
  } // end loop over throws
  
  if (par_throws) delete par_throws;
  
  wt_storer->SaveToFile(); // Save the weights to a file
  delete wt_storer;
  
#ifdef __T2KRW_OAANALYSIS_ENABLED__
  for(int i=0; i<(int)allRooTrackerVtxs.size(); ++i) delete allRooTrackerVtxs[i];
#endif
  
  for(int i=0; i<(int)allSKVtxs.size(); ++i) delete allSKVtxs[i];
  
#else
  
  cout << "The NEUT, JNuBeam and NIWG reweighting engines must all be enabled upon building." << endl;
  
#endif
  return 0;
}

// Print the cmd line syntax
void Usage(){
  cout << "Error Usage: " << endl;
  cout << "genWeights_2015.exe -i <Input Filename> -horn <+1:nu-mode, -1:antinu-mode> [-app <1 for appearance sample, 0 otherwise>] -p <Input BANFF Parameter Filename> -o <Output Weight Tree Filename>" << endl;
  exit(1);
}

// Messy way to process cmd line arguments.
void ParseArgs(int argc, char **argv){
  int nargs = 3;
  if(argc<(nargs*2+1)){ Usage(); }
  for(int i = 1; i < argc; i++){
    if(string(argv[i]) == "-i") {fInputFileName = argv[i+1]; i++;}
    else if(string(argv[i]) == "-horn") {iHornMode = std::atoi(argv[i+1]); i++;}
    else if(string(argv[i]) == "-app") {iApp = std::atoi(argv[i+1]); i++;}
    else if(string(argv[i]) == "-dslist") {fDisablSystList = argv[i+1]; i++;}
    else if(string(argv[i]) == "-o") {fOutputFileName = argv[i+1]; i++;}
    else if(string(argv[i]) == "-p") {fInputParsFileName = argv[i+1]; i++;}
    else if(string(argv[i]) == "-t") {nThrows = std::atoi(argv[i+1]); i++;}
    else if(string(argv[i]) == "-r") {rndSeed = std::atoi(argv[i+1]); i++;}
    else if(string(argv[i]) == "-e") {fNEvts = std::atoi(argv[i+1]); i++;}
    else if(string(argv[i]) == "--use-prefit") fUsePrefit = true;
    else if(string(argv[i]) == "--drop-flux") fUseFlux = false;
    else if(string(argv[i]) == "--drop-xsec") fUseXsec = false;
    else {
      cout << "Invalid argument:" << argv[i] << " "<< argv[i+1] << endl;
      Usage();
    }
  }
}

