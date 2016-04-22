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

#include "SK__h1.h"

using namespace std;
using namespace t2krew;

int fNskEvts = -1;
TString fSKFileName;
TString outfilename = "";
char *outhistofilename = "";

//added by maggie
int parset_defined=0;
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
    outfilename.ReplaceAll(".root","_weights.root");
  }
  
  // Parse command line input of systematic names and values
  vector<string> handles = separate(handles_input);
  vector<string> values = separate(values_input);
  
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
  
// TFile *outhistsfile = new TFile(outhistofilename,"RECREATE");

  /*
 const int nSamples = 2;
 TH1F *h_fsiwts[nSamples][2];
 TH2F *h_fsiwts_vs_nvert[nSamples][2];
 TH2F *h_fsiwts_vs_nvcvert[nSamples][2];

 TH1F *h_mode[nSamples];

 TString sampleName[nSamples] = {"MGeV 1R #nu_{e}", "MGeV MR #mu-like"};
 for (int isample=0; isample<nSamples; isample++) {
   
   h_fsiwts[isample][0] = new TH1F(Form("h_fsiwts_%d_min",isample),Form("%s; Min(FSI Weight)",sampleName[isample].Data()),150,-1,2);
   h_fsiwts[isample][1] = new TH1F(Form("h_fsiwts_%d_max",isample),Form("%s; Max(FSI Weight)",sampleName[isample].Data()),5000,0,200);

   h_fsiwts_vs_nvert[isample][0] = new TH2F(Form("h_fsiwts_vs_nvert_%d_min",isample),Form("%s; Min(FSI Weight)",sampleName[isample].Data()),150,-1,2,100,0,100);
   h_fsiwts_vs_nvert[isample][1] = new TH2F(Form("h_fsiwts_vs_nvert_%d_max",isample),Form("%s; Max(FSI Weight)",sampleName[isample].Data()),5000,0,200,100,0,100);

   h_fsiwts_vs_nvcvert[isample][0] = new TH2F(Form("h_fsiwts_vs_nvcvert_%d_min",isample),Form("%s; Min(FSI Weight)",sampleName[isample].Data()),150,-1,2,300,0,300);
   h_fsiwts_vs_nvcvert[isample][1] = new TH2F(Form("h_fsiwts_vs_nvcvert_%d_max",isample),Form("%s; Max(FSI Weight)",sampleName[isample].Data()),5000,0,200,300,0,300);
   
   h_mode[isample] = new TH1F(Form("h_mode_%d",isample),Form("%s;NEUT mode",sampleName[isample].Data()),106,-53,53);
   
 }
   */
  
  // Loop over parameter sets
  for (int iparset=0; iparset<nFSIsets; iparset++) {
    
    cout << "Var. #" << iFSI[iparset] << " : \t";
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
    
    // Set value for each FSI parameter
    rw->Systematics().SetTwkDial(t2krew::kNCasc_FrAbs_pi,     dial[FSIabs]);
    rw->Systematics().SetTwkDial(t2krew::kNCasc_FrCExLow_pi,  dial[FSIcxl]);
    rw->Systematics().SetTwkDial(t2krew::kNCasc_FrInelLow_pi, dial[FSIqel]);
    rw->Systematics().SetTwkDial(t2krew::kNCasc_FrPiProd_pi,  dial[FSIhad]);
    rw->Systematics().SetTwkDial(t2krew::kNCasc_FrCExHigh_pi, dial[FSIcxh]);
    rw->Systematics().SetTwkDial(t2krew::kNCasc_FrInelHigh_pi,dial[FSIqeh]);
    
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
      
      //if (!(skVtxs->evis > 30.0 && skVtxs->wall > 200.0 && skVtxs->nhitac < 16)) continue;
      //
      //int sampleType = -1;
      //
      //int nDecayE = muedcyBuild(skVtxs);
      //
      //// M.ost E.energetic R.ing index
      //int   mer;
      //float p;
      //float pmax  = 0.;
      //for( int i = 0 ; i < skVtxs->nring ; i++ ) {
      //  p = ( skVtxs->prmslg[i][1] < skVtxs->prmslg[i][2] ? skVtxs->amome[i] : skVtxs->amomm[i] );
      //  if( p > pmax ){ pmax = p; mer = i ;}
      //}
      //
      //// What is the identity of the M.ost E.nergetic R.ing
      //int mer_id = ( skVtxs->prmslg[mer][1] < skVtxs->prmslg[mer][2] ? 1 : 2 ); // 1:e 2:mu
      //
      //if (skVtxs->ip[0] != 3 && skVtxs->nring == 1 && skVtxs->evis >= 1330 && nDecayE>0) {
      //  sampleType = MultiGeV_elike_nue;
      //}
      //
      //else if (mer_id == 2 && skVtxs->nring > 1) {
      //  if (skVtxs->evis < 1330.0 && skVtxs->evis > 600.0 && pmax > 600.0 || skVtxs->evis > 1330.0) {
      //	 sampleType = MultiRing_mulike;
      //  }
      //}
      //
      //else  continue;
      //
      //if (sampleType<0) continue;
      //
      //
      //double fsi_wt_min =1;
      //double fsi_wt_max =1;
      //
      //for (int ipar=0; ipar<24; ipar++) {
      //  if(skVtxs->Fsivarwt[ipar] < fsi_wt_min) fsi_wt_min = skVtxs->Fsivarwt[ipar];
      //
      //  if(skVtxs->Fsivarwt[ipar] > fsi_wt_max) fsi_wt_max = skVtxs->Fsivarwt[ipar];
      //}
      //
      //
      //h_fsiwts[sampleType][0]->Fill(TMath::Max(-1.,fsi_wt_min));
      //h_fsiwts[sampleType][1]->Fill(TMath::Min(fsi_wt_max,199.99));
      //
      //h_fsiwts_vs_nvert[sampleType][0]->Fill(TMath::Max(-1.,fsi_wt_min),skVtxs->Nvert);
      //h_fsiwts_vs_nvert[sampleType][1]->Fill(TMath::Min(fsi_wt_max,199.99),skVtxs->Nvert);
      //
      //h_fsiwts_vs_nvcvert[sampleType][0]->Fill(TMath::Max(-1.,fsi_wt_min),skVtxs->Nvcvert);
      //h_fsiwts_vs_nvcvert[sampleType][1]->Fill(TMath::Min(fsi_wt_max,199.99),skVtxs->Nvcvert);
      //
      //if (!skVtxs->Nvert && skVtxs->Ibound==1) h_mode[sampleType]->Fill(skVtxs->mode);
      
      double weight = rw->CalcWeight(skVtxs);
      //if ((weight-skVtxs->Fsivarwt[iparset-1])/weight > 0.000001)
      //  cout << weight << " " << skVtxs->Fsivarwt[iparset-1] << endl;
      
      // Store weight in tree
      wt_storer->AddWeight(weight);
      
    }
    
  }
  
// outhistsfile->Write();

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

    case 'n':
      fNskEvts = atoi(argv[2]);
      ++argv; --argc;
      break;


    case 'o':
      outfilename = argv[2];
      ++argv; --argc;
      break;
      

    case 'h':
      outhistofilename = argv[2];
      ++argv; --argc;
      break;

    case 'v':
      values_input = argv[2];
      ++argv; --argc;
      break;

//added by maggie      
    case 'm':
      parset_defined = atoi(argv[2]);
      ++argv; --argc;
      break;
   
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


/*
 int muedcyBuild(SK::SK__h1 *skVtxs)
 {

   int skgen=3; // SK4

   int i , nmuemax;
   int lmuedcy = 0;

   //int ehit_cut_1[4] = { 60 , 30 , 60 , 60 };
   //int ehit_cut_2[4] = { 40 , 20 , 40 , 40 };

   nmuemax = ( skVtxs->nmue > 10 ? 10 :  skVtxs->nmue );
   if ( skVtxs->nmue < 0) return 0;

   for( i = 0 ; i < nmuemax ; i++ )
     {
       if(skgen ==3) //for SK4 only
	 {
	   if(skVtxs->etime[i]<0.1) continue;
	   lmuedcy++;
	 }
       //else
       //	 {
       //	   if(  skVtxs->evis <  1330.0 && skVtxs->etime(i) < 0.1 ) continue;
       //	   if(  skVtxs->evis >= 1330.0 && skVtxs->etime(i) < 1.2 ) continue;
       //	   if(  skVtxs->etime(i) > 0.8    && skVtxs->etime(i) < 1.2 ) continue;
       //
       //	   if( skVtxs->etype(i) == 1 && skVtxs->ehit(i) >= skVtxs->ehit_cut_1[skgen] && egood(i) > 0.5 )
       //	     {
       //	       lmuedcy++;
       //	     }
       //	   else if( skVtxs->etype(i) >= 2 && skVtxs->etype(i) <= 4 && skVtxs->ehit(i) >= skVtxs->ehit_cut_2[skgen] )
       //	     {
       //	       lmuedcy++;
       //	     }
       //	 }
     }// end of loop on dcy-e
   

   return lmuedcy;

}

*/
