//
// Example to show how to reweight NEUT events from the SK tree. Run with:
//
// ./genWeightsFromSK_NEUTFSI.exe -s "sk_inputfile(s)"
//                                -o output filename (without .root extension)
//                                -n <number of events to reweight>
//                                -h "list of systematic names"
//                                -v "corresponding list of systematic values"
//           //added by maggie
//                                -m define which set of parameters is used (0-24)
//                                -t tune multiplicity to CHORUS data
//
// Patrick de Perio - Nov. 4, 2011
//         This example has the FSI parameters enabled in the code (you do
//         not need to specify the paramters by command line, i.e. -h & -v).
//         It reads the FSI parameter values from a table (fsi_pars.txt) and converts
//         it to the proper units for T2KReWeight.
//
//
#include <stdlib.h>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"

#include "T2KReWeight.h"
#include "T2KSyst.h"
#include "T2KWeightsStorer.h"

#include "T2KNeutReWeight.h"
//#include "T2KNIWGReWeight.h"
//#include "T2KSKReWeight.h"
#include "T2KGEANTReWeight.h"

#include "SK__h1.h"

using namespace std;
using namespace t2krew;

int fNskEvts = -1;
TString fSKFileName;
TString outfilename = "";
//char *outhistofilename = "";

int flgRW = 0;

//added by maggie
//int parset_defined=0;
int chorus_tune_flag=0;

//enum sampleEnum {MultiGeV_elike_nue, MultiRing_mulike};

//int muedcyBuild(SK::SK__h1 *skVtxs);
//int mringtype(SK::SK__h1 *skVtxs);

int getArgs(int argc, char* argv[]);
vector<string> separate(string input_uptions);
void include_systematics(t2krew::T2KReWeight *rw, vector<string> handles, vector<string> values);
string handles_input="";
string values_input="";

int main(int argc, char *argv[])
{
  
  // process the arguments
  int args = getArgs(argc, argv);
  if(args != 0){
    cout << "Usage " << endl;
    return 0;
  }
  
  if (outfilename.CompareTo("")==0) {
    outfilename = fSKFileName;
    outfilename.ReplaceAll(".root","_wgt.root");
  }
  
  // Parse command line input of systematic names and values
  vector<string> handles = separate(handles_input);
  vector<string> values = separate(values_input);
  
  bool fDoFSIrw = true;
  bool fDoSIrw = true;
  if (flgRW%2) fDoSIrw = false;
  if ((flgRW/2)%2) fDoFSIrw = false;
  cout << "FSI reweight : " << ( fDoFSIrw ? "Yes" : "No" ) << endl;
  cout << "SI reweight : " << ( fDoSIrw ? "Yes" : "No" ) << endl;
  
  cout << "Starting to reweight NEUT events from SK file(s): " << fSKFileName << endl;
  
  
  // Load the tree(s)
  TChain *SK_tree = new TChain();
  SK_tree->Add(Form("%s/h1",fSKFileName.Data()));
  
  // Instantiate the tree object to be read by T2KReWeight
  SK::SK__h1 *skVtxs = new SK::SK__h1(SK_tree,1);
  
  // Make a t2kreweighting object and add a NEUT weighting engine.
  t2krew::T2KReWeight *rw = new T2KReWeight();
  rw->AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());
  //rw->AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
  //rw->AdoptWghtEngine("sk_rw", new t2krew::T2KSKReWeight());
  rw->AdoptWghtEngine("geant_rw", new t2krew::T2KGEANTReWeight());
  
  // For storing the tree of weights for each event
  T2KWeightsStorer *wt_storer = new T2KWeightsStorer(outfilename.Data());
  
  // Enable systematics from command line argument
  include_systematics(rw,handles,values);
  
  // Enable CHORUS tuning parameter
  //rw->Systematics().Include(t2krew::kNIWGPiMult_CorrSwitch);
  //rw->Systematics().SetAbsTwk(t2krew::kNIWGPiMult_CorrSwitch);
  //rw->Systematics().SetTwkDial(t2krew::kNIWGPiMult_CorrSwitch, chorus_tune_flag);
  
  // Enable NEUT FSI parameters
  rw->Systematics().Include(t2krew::kNCasc_FrAbs_pi);      // abs
  rw->Systematics().Include(t2krew::kNCasc_FrCExLow_pi);   // cxl
  rw->Systematics().Include(t2krew::kNCasc_FrInelLow_pi);  // qel
  rw->Systematics().Include(t2krew::kNCasc_FrPiProd_pi);   // had
  rw->Systematics().Include(t2krew::kNCasc_FrCExHigh_pi);  // cxh
  rw->Systematics().Include(t2krew::kNCasc_FrInelHigh_pi); // qeh
  
  // Set units to absolute fractional change from default value
  rw->Systematics().SetAbsTwk(t2krew::kNCasc_FrAbs_pi);
  rw->Systematics().SetAbsTwk(t2krew::kNCasc_FrCExLow_pi);
  rw->Systematics().SetAbsTwk(t2krew::kNCasc_FrInelLow_pi);
  rw->Systematics().SetAbsTwk(t2krew::kNCasc_FrPiProd_pi);
  rw->Systematics().SetAbsTwk(t2krew::kNCasc_FrCExHigh_pi);
  rw->Systematics().SetAbsTwk(t2krew::kNCasc_FrInelHigh_pi);
  
  rw->Systematics().Include(t2krew::kGEANT_PionXSecTbl);
  
  const int nFSIsets = 25;  // Number of FSI parameter variations
  const int nFSIpars = 6;   // Number of FSI parameters
  
  // Following enumeration must correspond to order of columns of input
  // FSI parameter value table (fsi_pars.txt; from Table 13, T2K-TN-032)
  enum FSIenum {FSIqel, FSIqeh, FSIhad, FSIabs, FSIcxl, FSIcxh};
  
  int iFSI[nFSIsets];
  double FSIval[nFSIsets][nFSIpars];  // Value of each FSI parameter
  
  // Load parameter values from file
  cout << "Opening loading FSI parameter variations..." << endl;
  ifstream fsipars("pi_fsi_tbl.txt");
  
  string strtemp;
  if( fsipars.is_open() ){
    
    getline(fsipars,strtemp);
    
    int iparset = 0;
    while ( !fsipars.eof() ) {
      
      if (iparset>=nFSIsets) {
        cout << "Too many rows in file!" << endl;
        break;
      }
      
      fsipars >> iFSI[iparset];
      for(int ipar=0; ipar<nFSIpars; ipar++){
        fsipars >> FSIval[iparset][ipar];
      }
      
      iparset++;
      
    }
    fsipars.close();
    
  }
  else {
    printf( "pi_fsi_tbl.txt file not found!!\n" );
    exit(-1);
  }
  
  cout << " Loaded " << nFSIsets << " variations:" << endl;
  cout << strtemp << endl;
  for(int iparset=0; iparset<nFSIsets; iparset++){
    cout << iFSI[iparset] << "\t";
    for(int ipar=0; ipar<nFSIpars; ipar++){
      cout << FSIval[iparset][ipar] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  
  
  // Loop over parameter sets
  for (int iparset=0; iparset<nFSIsets; iparset++) {
    
    cout << endl << "Var. #" << iFSI[iparset] << " : \t";
    for(int ipar=0; ipar<nFSIpars; ipar++){
      cout << FSIval[iparset][ipar] << "\t";
    }
    cout << endl;
    
//    if(iparset != parset_defined) continue;
    
    double dial[nFSIpars];
    // Convert FSI value to fractional change from default
    // (Warning: Assumes default is the first row of input table)
    for (int ipar=0; ipar<nFSIpars; ipar++) {
      dial[ipar] = (FSIval[iparset][ipar]-FSIval[0][ipar])/FSIval[0][ipar];
    }
    
    Double_t dialtbl = 0.; // for secondary int.
    if(iFSI[iparset]!=0) dialtbl = (iFSI[iparset]-7)%8+1;
    
    // Set value for each FSI parameter
    if (fDoFSIrw) {
      rw->Systematics().SetTwkDial(t2krew::kNCasc_FrAbs_pi,     dial[FSIabs]);
      rw->Systematics().SetTwkDial(t2krew::kNCasc_FrCExLow_pi,  dial[FSIcxl]);
      rw->Systematics().SetTwkDial(t2krew::kNCasc_FrInelLow_pi, dial[FSIqel]);
      rw->Systematics().SetTwkDial(t2krew::kNCasc_FrPiProd_pi,  dial[FSIhad]);
      rw->Systematics().SetTwkDial(t2krew::kNCasc_FrCExHigh_pi, dial[FSIcxh]);
      rw->Systematics().SetTwkDial(t2krew::kNCasc_FrInelHigh_pi,dial[FSIqeh]);
    }
    
    if (fDoSIrw) {
      rw->Systematics().SetTwkDial(t2krew::kGEANT_PionXSecTbl,dialtbl);
    }
    
    // Propagate changes to NEUT reweighting engine
    rw->Reconfigure();
    
    // Store the values of the current set of parameters
    wt_storer->NewSystSet(rw->Systematics());
    
    if(fNskEvts < 0) fNskEvts = skVtxs->GetEntries();
    
    // Loop over SK entries and calculate weights.
    cout << "Will reweight SK nevents: " << fNskEvts << endl;
    for(int i = 0; i < fNskEvts; i++){
      
      skVtxs->GetEntry(i);
      
      if (i%1000==0)
        cout << "Event " << i << endl;
      
      
      double weight = rw->CalcWeight(skVtxs);
      //if ((weight-skVtxs->Fsivarwt[iparset-1])/weight > 0.000001)
      //  cout << weight << " " << skVtxs->Fsivarwt[iparset-1] << endl;
      
      // Store weight in tree
      wt_storer->AddWeight(weight);
      
    }
    cout << " - Reweighting is complete!" << endl << endl;
    
  }
  
  wt_storer->SaveToFile();
  delete wt_storer;
  
  return 0;
}


int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){
    switch(argv[1][1]){

      // Get ROOT input root files
    case 's': 
      fSKFileName = argv[2];
      ++argv; --argc;
      break;

    case 'o':
      outfilename = argv[2];
      ++argv; --argc;
      break;
      
//    case 'h':
//      outhistofilename = argv[2];
//      ++argv; --argc;
//      break;

    case 'v':
      values_input = argv[2];
      ++argv; --argc;
      break;

    case 'r':
      flgRW = atoi(argv[2]);
      ++argv; --argc;
      break;

//added by maggie      
//    case 'm':
//      parset_defined = atoi(argv[2]);
//      ++argv; --argc;
//      break;
   
    case 't':
      chorus_tune_flag = atoi(argv[2]);
      ++argv; --argc;
      break;
    }

    ++argv; --argc;
  }

  return 0;

}


vector<string> separate(string input_options){
  vector<string> options_string;
  istringstream iss(input_options);
  while(iss){
    string sub;
    iss >> sub;
    options_string.push_back(sub);
  }
  options_string.pop_back();
  return options_string;
}

void include_systematics(t2krew::T2KReWeight *rw, vector<string>handles, vector<string>values) {
  
  // For reconfiguring the NEUT engine
  
  for(unsigned int i=0;i<handles.size();i++){
    rw->Systematics().Include(T2KSyst::FromString(handles[i])); 
    rw->Systematics().SetTwkDial(T2KSyst::FromString(handles[i]),atof(values[i].c_str()));
    
  }
  return;
}
