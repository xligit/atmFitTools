#ifndef visRing_C
#define visRing_C

#include "visRing.h"

void visRing::countdecaypi0(){
  return;
}

void visRing::countsecondary(){

  // rings to find in secondaries:
  // - gamma from pi0 decay
  // - photonuke or other non-brem gamma
  // - decay in pions

  // fqEvent has already been filled with current event info
  int ipid;
  double beta;
  nvisscnd = 0;
  
  for (int i=0; i<fq->nscndprt; i++){

     // decays near vertex
     if (fq->tscnd[i] > Tthresh) continue;

     // is decay and not scatter, etc
     if (fq->lmecscnd[i]!=5) continue;

     ipid=(int)TMath::Abs(fq->iprtscnd[i]); //< get ID code of particle
     
     // accept only certain decays
     int parentpid = TMath::Abs(fq->iorgprt[i]);
     if (ipid==parentpid) continue;
     if ((parentpid!=111) &&
        (parentpid!=211) &&
        (parentpid!=130) &&
        (parentpid!=310) &&
        (parentpid!=311) &&
        (parentpid!=3122) &&
        (parentpid!=13)) continue;

     // convert to geant code :(
     if (ipid==11) ipid=3;
     else if (ipid==13) ipid=5;
     else if (ipid==22) ipid=1;
     else if (ipid==211) ipid=9;
     else if (ipid==111) ipid=7;
     else{
       continue;
     }        
    
     // calc momentum
     double pmom = TMath::Sqrt(fq->pscnd[i][0]*fq->pscnd[i][0]+
                               fq->pscnd[i][1]*fq->pscnd[i][1]+
                               fq->pscnd[i][2]*fq->pscnd[i][2]);

     // see if above threshold
     if (ipid==1){
       // only count from pi0 decay
       if (TMath::Abs(fq->iorgprt[i])!=111) continue;
       if (pmom > gamthresh){
         addvisiblesecondary(ipid, i, pmom);
       }
     }
     else{
       // count
       beta = getbeta(ipid, pmom);
       if (beta<Cthresh) continue;
       else{
         addvisiblesecondary(ipid, i, pmom);
       }
     }

  }

  /*
  // loop over secondary particle stack
  for (int i=0; i<fq->nscndprt; i++){
      
    // make sure time is close to initial event
    if (fq->tscnd[i] > Tthresh) continue;

    ipid=(int)TMath::Abs(fq->iprtscnd[i]); //< get ID code of particle
    // convert to Geant code :(
    if (ipid==11) ipid=3;
    else if (ipid==13) ipid=5;
    else if (ipid==22) ipid=1;
    else if (ipid==211) ipid=7;
    else if (ipid==221) ipid=8;
    else if (ipid==2212) ipid=14;
    else{
      continue;
    }

    // for gamma
    if (ipid==1){ //< if particle is gamma, see if it will shower
      if ((fq->pmomv[i]>gamthresh)
      &&(fq->lmecscnd[i]<6)
      &&(fq->lmecscnd[i]>9)){ //< count as visible gamma if not from shower
        double pmom = TMath::Sqrt(fq->pscnd[i][0]*fq->pscnd[i][0]+
                                  fq->pscnd[i][1]*fq->pscnd[i][1]+
                                  fq->pscnd[i][2]*fq->pscnd[i][2]);
        addvisible(ipid, i, pmom);
      }
      else{
        continue; //< do nothing
      }
    }
    // for non-gaamma
    else{ 
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
      }
    }
  }
  */
  return;
}

double visRing::getpcrit(int ipid){
 double pcrit = 1e9;
 // use gamma threshold for gammas or neutral pions
 if ((ipid==1)||(ipid==7)){
   pcrit = gamthresh;
 }
 else{
   double mass = massof[ipid];
   pcrit =  TMath::Sqrt((Cthresh*Cthresh*mass*mass)/(1.-(Cthresh*Cthresh)));
 }
 return pcrit;
}

void visRing::addvisiblesecondary(int ipid, int index, double momentum){
  visscndpid[nvisscnd]=ipid;
  visscndparentid[nvisscnd] = fq->iorgprt[index];
  nvisscnd++;
  addvisible(ipid, index, momentum);
  return;
}

void visRing::addvisible(int ipid, int index, double momentum){

 
  // get critical momentum
  double mass = massof[ipid];
  double pcrit = getpcrit(ipid);

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
  if (ipid==7){
    pi0mom[nvpi0] = momentum;
    nvpi0++; //< pi0
   // visflg = 1;
  }

  if (visflg){
    // fill total arrays
    visindx[nvis]=index;
    vispid[nvis]=ipid;
    // for neutral pions, look at energy bound of the gamma
    if (ipid==7){
      double Epi0 = TMath::Sqrt(momentum*momentum + mass*mass);
      visstr[nvis] = ((Epi0/2.)-pcrit);
    }
    visstr[nvis] = momentum-pcrit;
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
  for (int i=2;i<fq->npar;i++){ //< outgoing particle index starts at 2
    ipid=(int)fq->ipv[i]; //< get ID code of particle
    if (ipid==1){ //< if particle is gamma, see if it will shower
      if ((fq->pmomv[i]>gamthresh)){
        addvisible(ipid, i, fq->pmomv[i]);
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
      }
    }
  }
  return;
}

void visRing::calcderived(){

  // find the maximum visble strength
  double maxstrength = 0.;
  int    imax = 0;
  for (int i=0; i<nvis; i++){
    if (visstr[i]>maxstrength){ maxstrength = visstr[i]; imax = i; }   
  }
  vismrpid1 = vispid[imax];

  // find second max visible strength
  double diff = 100000.;
  int iscnd = 0;
  for (int j=0; j<nvis; j++){
    if (j==imax) continue;
    double jdiff = TMath::Abs(visstr[j]-maxstrength);
    if (jdiff<diff) {iscnd = j; diff = jdiff; }
  }
  if (nvis==1) { vismrpar = 0.; vismrpid1=0; vismrpid2=0; }
  else if (nvis==0) { vismrpar = 0.; vismrpid1=0; vismrpid2=0;}
  else{
    vismrpar = visstr[iscnd];
    vismrpid2 = vispid[iscnd];
  }
  
  //
  return;
}

void visRing::testevent(int iev){
  cout<<"-----------------------------"<<endl;
  cout<<"--Primary--"<<endl;
  cout<<"N particles: "<<endl;
  for (int i=2; i<fq->npar; i++){
    double ipid = TMath::Abs((int)fq->ipv[i]);
    cout<<"  "<<ipid<<endl;
  }

  cout<<"--Visible--"<<endl;
  for (int i=0; i<nvis; i++){
    cout<<"  "<<vispid[i]<<endl;
  }

}

void visRing::fillVisVar(){
  
  // count rings in primary stack
  countprimary();



  // count rings in secondary stack
  countsecondary();


    // debuggin
   /*
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
  */

  // calculate derived quantities
  calcderived();

  return;
}

double visRing::getbeta(int ipid, double pmom){
  int pidcode = ipid;
  if (pidcode==0) return 0; //other thing?
  if (pidcode==4) return 0; //neutrino 
  if (pidcode==14) return 0; //neutron 
  if (pidcode==7) return 1; //pi0 always visible
  double mq = massof[pidcode]; // mass times absolute value of charge
  //cout<<"mass: "<<m<<endl;
  double E = sqrt(mq*mq + pmom*pmom);
  //cout<<"beta: "<<pmom/E<<endl;
  return (pmom/E);
}


#ifndef T2K
visRing::visRing(fqEvent* fqin){
  fq = fqin;

  // define some constants
  Cthresh = 0.7505; //Cherenkov threshold in c
  Tthresh = 130.; //cutoff time in ns to be counted in this event
  gamthresh = 10.; // min energy of gamma

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

  // define some constants
  Cthresh = 0.7505; //Cherenkov threshold in c
  Tthresh = 100.; //cutoff time in ns to be counted in this event
  gamthresh = 10.; // min energy of gamma

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

  //////////////////showerthresh = 78.33;
  //gammathresh  = 2*0.58;
}
#endif


#endif
