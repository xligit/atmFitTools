#ifndef visRing_C
#define visRing_C

#include "visRing.h"

void visRing::countdecaypi0(){
  return;
}

int visRing::hasdschild(int vcindex){
  
  int haschild = 0;
  for (int i=0; i<fq->Npvc; i++){
     int parentidx = (fq->Iorgvc[i] - 1);    
     if (parentidx==vcindex) haschild=1;
  }

  return haschild;
}

void visRing::countsecondary(){

  // rings to find in secondaries:
  // - gamma from pi0 decay
  // - photonuke or other non-brem gamma
  // - decay topions

  // fqEvent has already been filled with current event info
  
  int ipid;
  double beta;
  nvisscnd = 0; //< keep track of secondaries found in event
  
  for (int i=0; i<fq->nscndprt; i++){

     // decays near vertex
     if (fq->tscnd[i] > Tthresh) continue;

     // is decay and not scatter, etc
     if (fq->lmecscnd[i]!=5) continue;

     ipid=(int)TMath::Abs(fq->iprtscnd[i]); //< get ID code of particle
     
     // accept certain decays
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
     ipid = pdg2geant(ipid);
     if (ipid<0) continue;
    
     // calc momentum
     double pmom = TMath::Sqrt(fq->pscnd[i][0]*fq->pscnd[i][0]+
                               fq->pscnd[i][1]*fq->pscnd[i][1]+
                               fq->pscnd[i][2]*fq->pscnd[i][2]);

     // see if above threshold
     if (ipid==1){
       if (pmom > gamthresh){
         addvisiblesecondary(ipid, i, pmom,0);
       }
     }
     else{
       beta = getbeta(ipid, pmom);
       if (beta<Cthresh) continue;
       else{
         addvisiblesecondary(ipid, i, pmom,0);
       }
     }
  }
  
  //
  return;
}

double visRing::getpcrit(int ipid){
 double pcrit = 1e9;
 // use gamma threshold for gammas or neutral pions
 if (ipid==1){
   pcrit = gamthresh;
   return pcrit;
 }
 if (ipid==7){
   pcrit = 0.;
   return pcrit;
 }
 else{
   double mass = massof[ipid];
   pcrit =  (Cthresh*mass)/TMath::Sqrt((1.-(Cthresh*Cthresh)));
 }
 return pcrit;
}

double visRing::getEcrit(int ipid){
 double Ecrit = 1e9;
 // use gamma threshold for gammas or neutral pions
 if (ipid==1){
   Ecrit = gamthresh;
   return Ecrit;
 }
 if (ipid==7){
   Ecrit = 0.;
   return Ecrit;
 }
 else{
   double mass = massof[ipid];
   Ecrit =  (mass)/TMath::Sqrt((1.-(Cthresh*Cthresh)));
 }
 return Ecrit;
}


void visRing::addvisiblesecondary(int ipid, int index, double momentum, int flgverb){
  visscndpid[nvisscnd]=ipid;
  visscndparentid[nvisscnd] = fq->iorgprt[index];
  nvisscnd++;
  addvisible(ipid, index, momentum,0,1);

  // print out if necessary
  if (flgverb){
   cout<<"---------------------"<<endl;   
   cout<<" index:       "<<index<<endl;
   cout<<" pid:         "<<ipid<<endl;
   cout<<" parent:      "<<fq->iprntprt[index]<<endl;
   cout<<" origin:      "<<fq->iorgprt[index]<<endl;
   cout<<" process:     "<<fq->lmecscnd[index]<<endl;
   cout<<" parentrk:    "<<fq->iprnttrk[index]<<endl;
   cout<<" parentidx:   "<<fq->iprntidx[index]<<endl;
   cout<<" interaction: "<<fq->iflgscnd[index]<<endl;
  }


  return;
}

void visRing::addvisible(int ipid, int index, double momentum, int flgverb, int flgscnd){

 
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
    if (!flgscnd) vistime[nvis]=0.0;
    else{
      vistime[nvis] = fq->tscnd[index];
    }
    visstr[nvis] = momentum-pcrit;
    vismom[nvis] = momentum;
    nvis++;
  }

  if (flgverb){
    cout<<"========================"<<endl;
    cout<<"  visible? "<<visflg<<endl;
    cout<<"  index:   "<<index<<endl;
    cout<<"  escape?  "<<fq->Ichvc[index]<<endl;
    cout<<"  pid:     "<<ipid<<endl;
    cout<<"  pidpdg:  "<<fq->Ipvc[index]<<endl;
    cout<<"  mom:     "<<momentum<<endl;
    cout<<"  pcrit:   "<<pcrit<<endl;
    cout<<"  orgidx:  "<<fq->Iorgvc[index]<<endl;
    cout<<"  orgpart: "<<fq->Ipvc[ fq->Iorgvc[index] -1 ]<<endl;
    cout<<"  intcode: "<<fq->Iflvc[index]<<endl;
    cout<<"  orgcode: "<<fq->Iflvc[ fq->Iorgvc[index] -1 ]<<endl;
  }

  //
  return;
}

int  visRing::pdg2geant(int ipart){
  int abspid = TMath::Abs(ipart);
  if (abspid==2212){ return 14; }                                                                  
  if (abspid==2112){ return 13; }                                                                  
  if (abspid==3122){ return 18; }                                                                  
  if (abspid==310){ return 16; }                                                                  
  if (abspid==321){ return 11; }                                                                  
  if (abspid==221){ return 17; }                                                                  
  if (abspid==211){ return 8; }                                                                  
  if (abspid==111){ return 7; }                                                                  
  if (abspid==130){ return 10; }                                                                  
  if (abspid==22){ return 1; }                                                                  
  if (abspid==11){ return 2; }                                                                  
  if (abspid==12){ return 4; }                                                                  
  if (abspid==13){ return 6; }                                                                  
  if (abspid==14){ return 4; };
  return -1;

}
    
void visRing::initconstants(){

  // define some constants
  Cthresh = 0.7505; //Cherenkov threshold in c
  Tthresh = 50.; //cutoff time in ns to be counted in this event
  gamthresh = 10.; // min energy of gamma

  // particle masses in MeV
  massof[1] = 0;  //gamma
  massof[2] = 0.511; //positron
  massof[3] = 0.511; //electron 
  massof[4] = 0; //neutrino
  massof[5] = 105.7; //mu+
  massof[6] = 105.7; //mu-
  massof[7] = 134.98; //pi0
  massof[8] = 139.6; //pi+
  massof[9] = 139.6; //pi-
  massof[10] = 497.6; //K long
  massof[11] = 493.68; //K+
  massof[12] = 493.68; //K-
  massof[13] = 939.6; //n
  massof[14] = 938.3; //p
  massof[15] = 938.3; //p-
  massof[17] = 547.8; //eta
  massof[18] = 1115.6; //lambda

  return;
}

void visRing::countprimaryvc(){
  
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
  for (int i=2;i<fq->Npvc;i++){ //< outgoing particle index starts at 2
   
 //   cout<<"check part: "<<i<<endl;
    // make sure particle escapes
    if (fq->Ichvc[i]!=1) continue;

    // make sure particle has no children (from scatters etc.)
    if (hasdschild(i)) continue;

    // get pid
    int ipidpdg=(int)fq->Ipvc[i]; //< get ID code of particle
    ipid = pdg2geant(ipidpdg); //< convert to geant 
    if (ipid<0){
      cout<<"No conversion found for particle: "<<ipidpdg<<endl;
      continue;
    }
//    cout<<"pid: "<<ipid<<endl;;

    // particle kinematics
    double mass = massof[ipid];
    double momentum = fq->Abspvc[i];
    double energy   = TMath::Sqrt(mass*mass + momentum*momentum);

    if (ipid==1){ //< if particle is gamma, see if it will shower
      if ((momentum>gamthresh)){
        addvisible(ipid, i, momentum,0);
      }
      else{
        continue; //< do nothing
      }
    }
    else{ //count other rings
      beta = getbeta(ipid,momentum);
      if (beta<Cthresh){
        continue;
      }
      else{
        addvisible(ipid, i, momentum,0);
      }
    }
  }

  return;

}

void visRing::countprimary(){
  
  // fqEvent has already been filled with current event info
  int ipid;;
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
  if (nvis<=1) { vismrpar = 0.; vismrpid1=0; vismrpid2=0; vismrt1=0; vismrt2=0; }
  else{
    vismrpar = visstr[iscnd];
    vismrpid2 = vispid[iscnd];
    vismrt1 = vistime[imax];
    vismrt2 = vistime[iscnd];
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
//  countprimary();

  // count rings in primaryc VCWORK stack
  countprimaryvc();


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
  double mm = massof[pidcode]; // mass 
  double E = sqrt(mm*mm + pmom*pmom);
  return (pmom/E);
}


#ifndef T2K
visRing::visRing(fqEvent* fqin){
  fq = fqin;
  initconstants();
}
#endif
#ifdef T2K
visRing::visRing(t2kfqEvent* fqin){
  fq = fqin;
  initconstants();
}
#endif


#endif
