#ifndef visRing_C
#define visRing_C

#include "visRing.h"

void visRing::countsecondary(){

  // fqEvent has already been filled with current event info
  int ipid;
  double beta;
  //double Cthresh = 0.7505; //Cherenkov threshold in c
//  double Tthresh = 50.; //cutoff time in ns to be counted in this event
//  double gamthresh = 10.; // min energy of gamma
  nvisscnd = 0;
  for (int i=0; i<fq->nscndprt; i++){
      
    // make sure time is close to initial event
    if (fq->tscnd[i] > Tthresh) continue;
    if (fq->iprntidx[i]>0) continue;
    ipid=(int)TMath::Abs(fq->iprtscnd[i]); //< get ID code of particle
    // convert to Geant code :(
    if (ipid==11) ipid=3;
    else if (ipid==13) ipid=5;
    else if (ipid==22) ipid=1;
    else if (ipid==211) ipid=7;
    else if (ipid==2212) ipid=14;
    else{
      continue;
    }
    if (ipid==1){ //< if particle is gamma, see if it will shower
      if ((fq->pmomv[i]>gamthresh)
      &&(fq->lmecscnd[i]<6)
      &&(fq->lmecscnd[i]>9)){ //< count as visible gamma if not from shower
        double pmom = TMath::Sqrt(fq->pscnd[i][0]*fq->pscnd[i][0]+
                                  fq->pscnd[i][1]*fq->pscnd[i][1]+
                                  fq->pscnd[i][2]*fq->pscnd[i][2]);
        addvisible(ipid, i, pmom);
/*        visindx[nvis]=i;
        vispid[nvis]=ipid;
        visstr[nvis] = (pmom/gamthresh);
        vismom[nvis] = pmom;
        nvis++;
        nvisscnd++;
        nvgam++;
        */
      }
      else{
        continue; //< do nothing
      }
    }
    else{ //count non-shower rings
    double pmom = TMath::Sqrt(fq->pscnd[i][0]*fq->pscnd[i][0]+
                              fq->pscnd[i][1]*fq->pscnd[i][1]+
                              fq->pscnd[i][2]*fq->pscnd[i][2]);
      beta = getbeta(ipid,pmom);
      if (beta<Cthresh){
        continue;
      }
      else{
        // reject certain interaction modes
        if (fq->lmecscnd[i]==6) continue; //< no pair prod 
        if (fq->lmecscnd[i]==7) continue; //< no compton 
        if (fq->lmecscnd[i]==8) continue; //< no photo-electric 
        if (fq->lmecscnd[i]==9) continue; //< no brem 
        if (fq->lmecscnd[i]==12) continue; //< no hadron int. 
        if (fq->lmecscnd[i]==13) continue; //< no hadron coh. 
        if (fq->lmecscnd[i]==20) continue; //< no hadron inelast. 
        if (fq->lmecscnd[i]==21) continue; //< no muon capture 
        if (fq->lmecscnd[i]==30) continue; //< no below thresh 

        addvisible(ipid, i, pmom);

        /*
        int countflg = 0; //< flag to count as visible particle
        // visible muon
        if ((ipid==5)||(ipid==6)){ nvmu++; countflg=1;  }
        // visible electron
        if ((ipid==2)||(ipid==3)){ nve++; countflg=1;  }
        // visible charged pion
        if ((ipid==9)||(ipid==8)){ nvpip++; countflg=1;  }
        // visible proton
        if ((ipid==14)||(ipid==15)){ nvp++; countflg=1;  }
        // contains pi0
        if (ipid==7) nvpi0++;
        // count particles
        if (countflg){
          double mass = massof[ipid];
          double pcrit = TMath::Sqrt((Cthresh*Cthresh*mass*mass)/(1.-(Cthresh*Cthresh)));
          visindx[nvis]=i;
          vispid[nvis]=ipid;
          visstr[nvis] = pmom/pcrit;
          vismom[nvis] = pmom;
          nvis++;
          nvisscnd++;
        }
        */
      }
    }
  }
  return;
}

void visRing::addvisible(int ipid, int index, double momentum){

 
  // get critical momentum
  double mass = massof[ipid];
  double pcrit;
  if (ipid==1){
    // for gamma
    pcrit = gammathresh;
  }
  else{
    pcrit = TMath::Sqrt((Cthresh*Cthresh*mass*mass)/(1.-(Cthresh*Cthresh)));
  }

  // fill arrays for individual particle types
  int visflg = 0;
  if (ipid==1){ //< gamma
    gammom[nvgam] = momentum;
    nvgam++;
    visflg = 1;
  }
  if ((ipid==2)||(ipid==3)){ //< electron 
    emom[nve] = momentum;
    nve++;
    visflg = 1;
  }
  if ((ipid==5)||(ipid==6)){ //< muon 
    mumom[nvmu] = momentum;
    nvmu++;
    visflg = 1;
  }
  if ((ipid==8)||(ipid==9)){ //< charged pion 
    pimom[nvpip] = momentum;
    nvpip++;
    visflg = 1;
  }
  if ((ipid==14)||(ipid==15)){ //< proton 
    protmom[nvp] = momentum;
    nvp++;
    visflg = 1;
  }
  if (ipid==7) nvpi0++; //< pi0

  if (visflg){
    // fill total arrays
    visindx[nvis]=index;
    vispid[nvis]=ipid;
    visstr[nvis] = (momentum/pcrit);
    vismom[nvis] = momentum;
    nvis++;
  }

  //
  return;
}

void visRing::countprimary(){
  
  // fqEvent has already been filled with current event info
  int ipid;
  double beta;

  // set initial visble particles to zero
  nvmu=0;
  nve=0;
  nvpip=0;
  nvpi0=0;
  nvoth=0;
  nvis=0;
  nvgam=0;
  nvp=0;
  nvk=0;
  
  // loop over primary particles
  for (int i=0;i<fq->npar;i++){
    ipid=(int)fq->ipv[i]; //< get ID code of particle
    if (ipid==1){ //< if particle is gamma, see if it will shower
      double gamthresh = 10.0;
      if ((fq->pmomv[i]>gamthresh)){
        addvisible(ipid, i, fq->pmomv[i]);
        // count as visible gamma 
    //    visindx[nvis]=i;
    //    vispid[nvis]=ipid;
    //    visstr[nvis] = (fq->pmomv[i]/gamthresh);
    //    vismom[nvis] = fq->pmomv[i];
    //    nvis++;
    //    nvgam++;
      }
      else{
        continue; //< do nothing
      }
    }
    else{ //count non-shower rings
      beta = getbeta(ipid,fq->pmomv[i]);
      if (beta<Cthresh){
        continue;
      }
      else{
        addvisible(ipid, i, fq->pmomv[i]);
        /*
        int countflg = 0; //< flag to count as visible particle
        // visible muon
        if ((ipid==5)||(ipid==6)){ nvmu++; countflg=1;  }
        // visible electron
        if ((ipid==2)||(ipid==3)){ nve++; countflg=1;  }
        // visible charged pion
        if ((ipid==9)||(ipid==8)){ nvpip++; countflg=1;  }
        // visible kaon
        //if ((ipid==11)||(ipid==12)){ nvk++; countflg=1;  }
        // visible proton
        if ((ipid==14)||(ipid==15)){ nvp++; countflg=1;  }
        // contains pi0
        if (ipid==7) nvpi0++;
        // count particles
        if (countflg){
          double mass = massof[ipid];
          double pcrit = TMath::Sqrt((Cthresh*Cthresh*mass*mass)/(1.-(Cthresh*Cthresh)));
          visindx[nvis]=i;
          vispid[nvis]=ipid;
          visstr[nvis] = fq->pmomv[i]/pcrit;
          if (visstr[nvis] > 10000.) cout<<"!!!" <<fq->pmomv[i]<<endl;
          vismom[nvis] = fq->pmomv[i];
          nvis++;
        }
        */
      }
    }
  }
  
  //
  return;
}

void visRing::calcderived(){

  // find the maximum visble strength
  double maxstrength = 0.;
  int    imax = 0;
  for (int i=0; i<nvis; i++){
    if (visstr[i]>maxstrength){ maxstrength = visstr[i]; imax = i; }   
  }
  // find second max visible strength
  double diff = 100000.;
  int iscnd = 0;
  for (int j=0; j<nvis; j++){
    if (j==imax) continue;
    double jdiff = TMath::Abs(visstr[j]-maxstrength);
    if (jdiff<diff) {iscnd = j; diff = jdiff; }
  }
  if (nvis==1) vismrstr = 1.;
  else if (nvis==0) vismrstr = 0.;
  else{
    vismrstr = 1. + visstr[iscnd]/maxstrength;
  }
  
  //
  return;
}

void visRing::fillVisVar(){
  
  // count rings in primary stack
  countprimary();

  // debuggin
  if (nvis==0) {
    for (int i=0; i<fq->npar; i++){
      int ipid = (int)fq->ipv[i];
      if ((ipid==5)||(ipid==6)){
        cout<<"mu mom: "<<fq->pmomv[i]<<endl;
      }
      if ((ipid==8)||(ipid==9)){
        cout<<"pi mom: "<<fq->pmomv[i]<<endl;
      }
      if ((ipid==2)||(ipid==3)){
        cout<<"e mom: "<<fq->pmomv[i]<<endl;
      }
      if (ipid==7){
        cout<<"pi0 mom: "<<fq->pmomv[i]<<endl;
      }
    }
  }

  // count rings in secondary stack
  countsecondary();

  // calculate derived quantities
  calcderived();

  return;
}

double visRing::getbeta(int ipid, double pmom){
  if (ipid==0) return 0; //other thing?
  if (ipid==4) return 0; //neutrino 
  if (ipid==14) return 0; //neutron 
  if (ipid==7) return 1; //pi0 always visible
  double mq = massof[ipid]; // mass times absolute value of charge
  //cout<<"mass: "<<m<<endl;
  double E = sqrt(mq*mq + pmom*pmom);
  //cout<<"beta: "<<pmom/E<<endl;
  return (pmom/E);
}

#ifndef T2K
visRing::visRing(fqEvent* fqin){
  fq = fqin;
  massof[1] = 0;  //gamma
  massof[2] = 0.511e1; //positron
  massof[3] = 0.511e1; //electron 
  massof[4] = 1e10; //neutrino
  massof[5] = 0.1057e3; //mu+
  massof[6] = 0.1057e3; //mu-
  massof[7] = 0.1350e3; //pi0
  massof[8] = 0.1396e3; //pi+
  massof[9] = 0.1396e3; //pi-
//  massof[10] = 0.4977e3; //K long
  massof[10] = 1e10; //K long
  massof[11] = 0.4937e3; //K+
  massof[12] = 0.4937e3; //K-
//  massof[13] = 0.9396e3; //n
  massof[13] = 1e10; //n
  massof[14] = 0.9383e3; //p
  massof[15] = 0.9383e3; //p-
  massof[17] = 0.5475e3; //eta
}
#endif
#ifdef T2K
visRing::visRing(t2kfqEvent* fqin){
  fq = fqin;
  massof[1] = 0;  //gamma
  massof[2] = 0.511e1; //positron
  massof[3] = 0.511e1; //electron 
  massof[4] = 1e10; //neutrino
  massof[5] = 0.1057e3; //mu+
  massof[6] = 0.1057e3; //mu-
  massof[7] = 0.1350e3; //pi0
  massof[8] = 0.1396e3; //pi+
  massof[9] = 0.1396e3; //pi-
//  massof[10] = 0.4977e3; //K long
  massof[10] = 1e10; //K long
  massof[11] = 0.4937e3; //K+
  massof[12] = 0.4937e3; //K-
//  massof[13] = 0.9396e3; //n
  massof[13] = 1e10; //n
  massof[14] = 0.9383e3; //p
  massof[15] = 0.9383e3; //p-
  massof[17] = 0.5475e3; //eta

  showerthresh = 78.33;
  gammathresh  = 2*0.58;
}
#endif

#endif
