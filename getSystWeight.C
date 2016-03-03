#ifndef GETSYSTWEIGHT_C
#define GETSYSTWEIGHT_C
#include "atmFitPars.C"
#include "fQreader.C"
#include "TString.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// A function to calculate the systematic error weight that should be assigned 
// to this event.
double getSystWeight(const char* systype, fQreader* mcevent, int ipar, double value){
  
  //initial weight is always one:
  double ww = 1.0;

  //get string to determine the parameterization type
  TString parameterization_type = systype;

  /////////////////////////////// 
  //simple TN186 parameterization
  if (!parameterization_type.CompareTo("tn186")){ 
    
    //these values are needed to determine the event weight
    int absmode = TMath::Abs(mcevent->mode);
    double Enu     = mcevent->pmomv[0];
    int  nutype  = TMath::Abs(mcevent->ipnu[0]);  

    //CCQE norm bin1 
    if (ipar==0){
      if ((absmode==1)&&(Enu<200.)) ww*=value;
    }
    //CCQE norm bin2 
    if (ipar==1){
      if ((absmode==1)&&(Enu>200.)&&(Enu<400.)) ww*=value;
    }
    //CCQE norm bin3 
    if (ipar==2){
      if ((absmode==1)&&(Enu>400.)&&(Enu<800.)) ww*=value;
    }
    //CCQE norm bin4 
    if (ipar==3){
      if ((absmode==1)&&(Enu>800.)) ww*=value;
    }
    //SubGevFlux
    if (ipar==4){
      if (Enu<1000.) ww*=value;
    }
    //MultiGeVFlux
    if (ipar==5){
      if (Enu>1000.) ww*=value;
    }
    //CCnQE
    if (ipar==6){
      if ((absmode>1)&&(absmode<30)) ww*=value;
    }
    //NC
    if (ipar==7){
      if (absmode>=30) ww*=value;
    }
    //mu2e ratio
    if (ipar==8){
      if (nutype==14) ww*=value;
    }
  }

  /////////////////////////////// 
  //cosmic muons systematics
  if (!parameterization_type.CompareTo("cosmic")){ 
    
    //these values are needed to determine the event weight
    int fvbin = mcevent->nbin; //< fiducial bin of event

    //FV Bin 1
    if (ipar==0){
      if (fvbin==0) ww*=value;
    }
    //FV Bin 2
    if (ipar==1){
      if (fvbin==1) ww*=value;
    }
    //FV Bin 3
    if (ipar==2){
      if (fvbin==2) ww*=value;
    }

  }

  //no negative weights
  if (ww<0.) ww = 0.;

  //////////
  return ww;


}



#endif
