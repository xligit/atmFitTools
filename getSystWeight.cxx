#ifndef GETSYSTWEIGHT_C
#define GETSYSTWEIGHT_C

#include "atmFitPars.cxx"
#include "fqProcessedEvent.h"

#include "TString.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// A function to calculate the systematic error weight that should be assigned 
// to this event give the value "value" of parameter "ipar"
double getSystWeight(const char* systype, fqProcessedEvent* mcevent, int ipar, double value){
  
  //initial weight is always one:
  double ww = 1.0;

  //get string to determine the parameterization type
  TString parameterization_type = systype;

  //simple TN186 parameterization
  if (!parameterization_type.CompareTo("tn186simple")){ 
    
    //these values are needed to determine the event weight
    int absmode = TMath::Abs(mcevent->mode); //< NEUT mode for this event
    double Enu     = mcevent->pmomv[0]; //< True neutrino Energy
    int  nutype  = TMath::Abs(mcevent->ipnu[0]); //< True neutrino type 

    // all of these parameters are sclaing parameters. If the value is <0,
    // then just set weight to 0 
    if (value<0.){
      ww=0.;
      return ww;
    }

    //CCQE norm bin1 
    if (ipar==0){
      if ((absmode==1)&&(Enu<190.)) ww*=value;
    }
    //CCQE norm bin2 
    if (ipar==1){
      if ((absmode==1)&&(Enu>190.)&&(Enu<240.)) ww*=value;
    }
    //CCQE norm bin3 
    if (ipar==2){
      if ((absmode==1)&&(Enu>240.)&&(Enu<294.)) ww*=value;
    }
    //CCQE norm bin4 
    if (ipar==3){
      if ((absmode==1)&&(Enu>294.)&&(Enu<333.)) ww*=value;
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



  //full TN186 parameterization
  if (!parameterization_type.CompareTo("tn186")){ 
    
    //these values are needed to determine the event weight
    int absmode = TMath::Abs(mcevent->mode); //< NEUT mode for this event
    double Enu     = mcevent->pmomv[0]; //< True neutrino Energy
    int  nutype  = TMath::Abs(mcevent->ipnu[0]); //< True neutrino type 

    // all of these parameters are sclaing parameters. If the value is <0,
    // then just set weight to 0 
    if (value<0.){
      ww=0.;
      return ww;
    }

    //CCQE norm bin1 
    if (ipar==0){
      if ((absmode==1)&&(Enu<190.)) ww*=value;
    }
    //CCQE norm bin2 
    if (ipar==1){
      if ((absmode==1)&&(Enu>190.)&&(Enu<240.)) ww*=value;
    }
    //CCQE norm bin3 
    if (ipar==2){
      if ((absmode==1)&&(Enu>240.)&&(Enu<294.)) ww*=value;
    }
    //CCQE norm bin4 
    if (ipar==3){
      if ((absmode==1)&&(Enu>294.)&&(Enu<333.)) ww*=value;
    }
    //CCQE norm bin5 
    if (ipar==4){
      if ((absmode==1)&&(Enu>333.)&&(Enu<390.)) ww*=value;
    }
    //CCQE norm bin6 
    if (ipar==5){
      if ((absmode==1)&&(Enu>390.)&&(Enu<440.)) ww*=value;
    }
    //CCQE norm bin7 
    if (ipar==6){
      if ((absmode==1)&&(Enu>440.)&&(Enu<487.)) ww*=value;
    }
    //CCQE norm bin8 
    if (ipar==7){
      if ((absmode==1)&&(Enu>487.)&&(Enu<590.)) ww*=value;
    }
    //CCQE norm bin9 
    if (ipar==8){
      if ((absmode==1)&&(Enu>590.)&&(Enu<690.)) ww*=value;
    }
    //CCQE norm bin10 
    if (ipar==9){
      if ((absmode==1)&&(Enu>690.)&&(Enu<786.)) ww*=value;
    }
    //CCQE norm bin11 
    if (ipar==10){
      if ((absmode==1)&&(Enu>786.)&&(Enu<896.)) ww*=value;
    }
    //CCQE norm bin12 
    if (ipar==11){
      if ((absmode==1)&&(Enu>896.)&&(Enu<994.)) ww*=value;
    }
    //CCQE norm bin13 
    if (ipar==12){
      if ((absmode==1)&&(Enu>994.)&&(Enu<2000.)) ww*=value;
    }
    //CCQE norm bin14 
    if (ipar==13){
      if ((absmode==1)&&(Enu>2000.)&&(Enu<3000.)) ww*=value;
    }
    //SubGevFlux
    if (ipar==14){
      if (Enu<1000.) ww*=value;
    }
    //MultiGeVFlux
    if (ipar==15){
      if (Enu>1000.) ww*=value;
    }
    //CCnQE
    if (ipar==16){
      if ((absmode>1)&&(absmode<30)) ww*=value;
    }
    //NC
    if (ipar==17){
      if (absmode>=30) ww*=value;
    }
    //mu2e ratio
    if (ipar==18){
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
