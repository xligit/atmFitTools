#ifndef visRing_C
#define visRing_C

#include "visRing.h"

void visRing::countdecaypi0(){
  return;
}


double visRing::getVisTowall(int index, int scndflg){

  TVector3 vpos;
  TVector3 vmom;
  TVector3 vdir;

  if (scndflg){
    vpos.SetXYZ(fq->vtxscnd[index][0],
                fq->vtxscnd[index][1],
                fq->vtxscnd[index][2]);
    vmom.SetXYZ(fq->pscnd[index][0],
                fq->pscnd[index][1],
                fq->pscnd[index][2]);
    vdir = vmom.Unit();
  }
  else{
    vpos.SetXYZ(fq->posv[0],
                fq->posv[1],
                fq->posv[2]);
    vmom.SetXYZ(fq->Pvc[index][0],
                fq->Pvc[index][1],
                fq->Pvc[index][2]);
    vdir = vmom.Unit();
  }
 
  double towallvalue = calcToWall(&vpos, &vdir);
  return towallvalue;
}

double visRing::getVisWall(int index, int scndflg){

  TVector3 vpos;
  TVector3 vmom;
  TVector3 vdir;

  if (scndflg){
    vpos.SetXYZ(fq->vtxscnd[index][0],
                fq->vtxscnd[index][1],
                fq->vtxscnd[index][2]);
    vmom.SetXYZ(fq->pscnd[index][0],
                fq->pscnd[index][1],
                fq->pscnd[index][2]);
    vdir = vmom.Unit();
  }
  else{
    vpos.SetXYZ(fq->posv[0],
                fq->posv[1],
                fq->posv[2]);
    vmom.SetXYZ(fq->Pvc[index][0],
                fq->Pvc[index][1],
                fq->Pvc[index][2]);
    vdir = vmom.Unit();
  }
 
  double wallvalue = calcWall(&vpos);
  return wallvalue;

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

     // reject certiain interactions
     if (fq->lmecscnd[i]==6) continue;  //< pair prod.
     if (fq->lmecscnd[i]==7) continue;  //< compton
     if (fq->lmecscnd[i]==8) continue;  //< photoelectric
     if (fq->lmecscnd[i]==9) continue;  //< brem
     if (fq->lmecscnd[i]==12) continue; //< hadronic int.
     if (fq->lmecscnd[i]==13) continue; //< coh. scatt.
     if (fq->lmecscnd[i]==30) continue; //< below thresh

     ipid=(int)TMath::Abs(fq->iprtscnd[i]); //< get ID code of particle
     
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
         addvisiblesecondary(ipid, i, pmom);
       }
     }
     else{
       beta = getbeta(ipid, pmom);
       if (beta<Cthresh) continue;
       else{
         addvisiblesecondary(ipid, i, pmom);
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
   Ecrit = gamthresh;
   return Ecrit;
 }
 else{
   double mass = massof[ipid];
   Ecrit =  (mass)/TMath::Sqrt((1.-(Cthresh*Cthresh)));
 }
 return Ecrit;
}

double visRing::getEnergy(int ipid, double pmom){
  double mm = massof[ipid];
  double E = TMath::Sqrt( mm*mm + pmom*pmom);
  return E;
}

double visRing::getKEcrit(int ipid){
   // use gamma threshold for gammas or neutral pions
 double Tcrit = 0.;
 if (ipid==1){
   Tcrit = gamthresh;
   return Tcrit;
 }
 if (ipid==7){
   Tcrit = 0;
   return Tcrit;
 }
 double mass = massof[ipid];
 Tcrit = (mass/TMath::Sqrt(1 - (Cthresh* Cthresh))) - mass;
 return Tcrit;
}

double visRing::getKE(int ipid, double pmom){
  if (ipid==1){
    return pmom;
  }
  double E = getEnergy(ipid, pmom);
  double m = massof[ipid];
  return E-m;
}

void visRing::addvisiblesecondary(int ipid, int index, double momentum){
  visscndpid[nvisscnd]=ipid;
  visscndparentid[nvisscnd] = fq->iorgprt[index];
  nvisscnd++;
  int visflg = addvisible(ipid, index, momentum,1);

  // print out if necessary
  if (flgverbscnd && visflg){
   cout<<"---------------------"<<endl;   
   int gpid = pdg2geant(fq->iprtscnd[index]);
   cout<<" pid:            "<<nameof[gpid]<<endl;
   cout<<" track:          "<<fq->itrkscnd[index]<<endl;
   cout<<" stack:          "<<fq->istakscnd[index]<<endl;
   cout<<" momentum:       "<<momentum<<endl;
   cout<<" time:           "<<fq->tscnd[index]<<endl;
   cout<<" interaction:    "<<fq->lmecscnd[index]<<endl;
   cout<<" parent id:      "<<fq->iprntprt[index]<<endl;
   cout<<" parent trackid: "<<fq->iprnttrk[index]<<endl;
   cout<<" parent index:   "<<fq->iprntidx[index]<<endl;
   cout<<" parent orig:    "<<fq->iorgprt[index]<<endl;
   cout<<" iflag:          "<<fq->iflgscnd[index]<<endl;
   cout<<" childs:         "<<fq->nchilds[index]<<endl;
 }


  return;
}

double visRing::getVisibleEnergy(int ipid, double pmom){

//  double E = getEnergy(ipid, pmom);

  //double Ecrit = getEcrit(ipid);

  double KE = getKE(ipid, pmom);
  double KEc = getKEcrit(ipid);
  return KE-KEc;
}

int visRing::addvisible(int ipid, int index, double momentum, int flgscnd){

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
    // fill arrays
    visindx[nvis]=index;
    vispid[nvis]=ipid;
    viswall[nvis] = getVisWall(index,flgscnd);
    vistowall[nvis] = getVisTowall(index,flgscnd);
    // for neutral pions, look at energy bound of the gamma
    if (ipid==7){
      visbrightness[nvis] = getVisibleEnergy(ipid,momentum);
    }
    if (!flgscnd) vistime[nvis]=0.0;
    else{
      vistime[nvis] = fq->tscnd[index];
    }
    visbrightness[nvis] = getVisibleEnergy(ipid,momentum);
    vismom[nvis] = momentum;
    nvis++;
  }

  if (flgverbprime){
    cout<<"========================"<<endl;
    cout<<"  visible? "<<visflg<<endl;
    cout<<"  index:   "<<index<<endl;
    cout<<"  escape?  "<<fq->Ichvc[index]<<endl;
    cout<<"  pid:     "<<ipid<<endl;
    cout<<"  pidpdg:  "<<fq->Ipvc[index]<<endl;
    cout<<"  mom:     "<<momentum<<endl;
    cout<<"  orgidx:  "<<fq->Iorgvc[index]<<endl;
    cout<<"  orgpart: "<<fq->Ipvc[ fq->Iorgvc[index] -1 ]<<endl;
    cout<<"  intcode: "<<fq->Iflvc[index]<<endl;
    cout<<"  orgcode: "<<fq->Iflvc[ fq->Iorgvc[index] -1 ]<<endl;
  }

  //
  return visflg;
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
  Tthresh = 10.; //cutoff time in ns to be counted in this event
  gamthresh = 10.; // min energy of gamma
  flgverbprime = 0;
  flgverbscnd= 0;
;
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

  // particle names (geant codes)
  nameof[1] = "gamma";
  nameof[2] = "positron";
  nameof[3] = "electron";
  nameof[4] = "neutrino";
  nameof[5] = "muon";
  nameof[6] = "antimuon";
  nameof[7] = "pi0";
  nameof[8] = "pi+";
  nameof[9] = "pi-";
  nameof[10] = "KL";
  nameof[11] = "K+";
  nameof[12] = "K-";
  nameof[13] = "n";
  nameof[14] = "p";
  nameof[15] = "p-";
  nameof[16] = "KS";
  nameof[17] = "eta";
  nameof[18] = "lambda";;
  
  
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
   
    // make sure particle escapes
    if (fq->Ichvc[i]!=1) continue;

    // make sure particle has no children (want most downstream particle);
    if (hasdschild(i)) continue;

    // get pid
    int ipidpdg=(int)fq->Ipvc[i]; //< get ID code of particle
    ipid = pdg2geant(ipidpdg); //< convert to geant 
    if (ipid<0){
      cout<<"No conversion found for particle: "<<ipidpdg<<endl;
      continue;
    }

    // particle kinematics
    double mass = massof[ipid];
    double momentum = fq->Abspvc[i];

    if (ipid==1){ //< if particle is gamma, see if it will shower
      if ((momentum>gamthresh)){
        addvisible(ipid, i, momentum);
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
        addvisible(ipid, i, momentum);
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

    if (fq->light_flag[i]!=1) continue;

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
    if (visbrightness[i]>maxstrength){ maxstrength = visbrightness[i]; imax = i; }   
  }
  vismrpid1 = vispid[imax];

  // find second max visible strength
  double diff = 100000.;
  int iscnd = 0;
  for (int j=0; j<nvis; j++){
    if (j==imax) continue;
    double jdiff = TMath::Abs(visbrightness[j]-maxstrength);
    if (jdiff<diff) {iscnd = j; diff = jdiff; }
  }
  if (nvis<=1) { vismrbrightness = 0.; vismrpid2=0; vismrt1=0; vismrt2=0; }
  else{
    vismrbrightness = visbrightness[iscnd];
    vismrpid2 = vispid[iscnd];
    vismrt1 = vistime[imax];
    vismrt2 = vistime[iscnd];
  }

  // get ring types
  if (vismrpid1==1 || vismrpid1==2 || vismrpid1==3){
    vismrtype1 = 1; 
  }
  else{
    vismrtype1 = 0;
  }
  if (vismrpid2==1 || vismrpid2==2 || vismrpid2==3){
    vismrtype2 = 1; 
  }
  else{
    vismrtype2 = 0;
  }

  // get FV info
  vismrwall1 = viswall[imax];
  vismrtowall1 = vistowall[imax];
  vismrwall2 = viswall[iscnd];
  vismrtowall2 = vistowall[iscnd];
  vismrwallmin = TMath::Min(vismrwall1, vismrwall2);
  vismrtowallmin = TMath::Min(vismrtowall1, vismrtowall2);
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

  // count rings in primaryc VCWORK stack
  countprimaryvc();

  // count rings in secondary stack
  countsecondary();

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

void visRing::printsecondaryindex(int ipid){
  for (int i=0; i<fq->nscndprt; i++){
    int abspid = TMath::Abs(fq->iprtscnd[i]);
    if (abspid==ipid){
      cout<<"found "<<ipid<<" at index: "<<i<<endl;
    }
  }
  return;
}


void visRing::printsecondaryinfo(int indx){
   cout<<"------------------------------"<<endl;
     int i =indx;
     double px = fq->pscnd[i][0];
     double py = fq->pscnd[i][1];
     double pz = fq->pscnd[i][2];
     double momentum = TMath::Sqrt(px*px + py*py + pz*pz);
     cout<<" pid:            "<<fq->iprtscnd[i]<<endl;
     cout<<" track:          "<<fq->itrkscnd[i]<<endl;
     cout<<" stack:          "<<fq->istakscnd[i]<<endl;
     cout<<" momentum:       "<<momentum<<endl;
     cout<<" time:           "<<fq->tscnd[i]<<endl;
     cout<<" interaction:    "<<fq->lmecscnd[i]<<endl;
     cout<<" parent id:      "<<fq->iprntprt[i]<<endl;
     cout<<" parent trackid: "<<fq->iprnttrk[i]<<endl;
     cout<<" parent index:   "<<fq->iprntidx[i]<<endl;
     cout<<" parent orig:    "<<fq->iorgprt[i]<<endl;
     cout<<" iflag:          "<<fq->iflgscnd[i]<<endl;
     cout<<" childs:         "<<fq->nchilds[i]<<endl;

   return; 
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
