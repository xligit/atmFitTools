#ifndef visRing_C
#define visRing_C

#include "visRing.h"

void visRing::fillVisVarFV(){
  //assume fqReader has already been filled with current event info
  int ipid;
  float beta;
  float thresh = 1./1.3333333;
  nvmu=0;
  nve=0;
  nvpip=0;
  nvpi0=0;
  nvoth=0;
  nvis=0;
  nvgam=0;
  nvk=0;
  nvp=0;
//  int nvisarr[MAXNVIS];
  for (int jpid=0;jpid<=50;jpid++){
    nvisarr[jpid]=0;
  }
  for (int i=0;i<fqfv->npar;i++){
    ipid=(int)fqfv->ipv[i];
    beta = getbeta(ipid,fqfv->pmomv[i]);
    if (beta<thresh) continue;
    else{
      visindx[nvis]=i;
      vispid[nvis]=ipid;
      nvis++;
      nvisarr[ipid]++;
      if ((ipid==5)||(ipid==6)) nvmu++;
      if ((ipid==2)||(ipid==3)) nve++;
      if ((ipid==9)||(ipid==8)) nvpip++;
      if ((ipid==11)||(ipid==12)) nvk++;
      if ((ipid==14)||(ipid==15)) nvp++;
      if ((ipid==7)) nvpi0++;
      if (ipid==1) nvgam++;
    }
  }
//  nvmu = nvisarr[5]+nvisarr[6];
//  nve  = nvisarr[2]+nvisarr[3];
//  nvpip = nvisarr[9]+nvisarr[8];
//  nvpi0 = nvisarr[7];
//  nvgam = nvisarr[1];
//  nvp   = nvisarr[14]+nvisarr[15];
//  nvk   = nvisarr[11]+nvisarr[12];
  nvoth=nvis-nvmu-nve-nvpip-nvpi0-nvgam-nvp-nvk;
  return;
}

void visRing::fillVisVar(){
  if (FVflg){
    fillVisVarFV();
    return;
  }
  //assume fqReader has already been filled with current event info
  int ipid;
  float beta;
  float thresh = 1./1.3333333;
  nvmu=0;
  nve=0;
  nvpip=0;
  nvpi0=0;
  nvoth=0;
  nvis=0;
  nvgam=0;
  nvp=0;
  nvk=0;
//  int nvisarr[MAXNVIS];
//  for (int jpid=0;jpid<=50;jpid++){
//    nvisarr[jpid]=0;
//  }
  for (int i=0;i<fq->npar;i++){
    ipid=(int)fq->ipv[i];
    beta = getbeta(ipid,fq->pmomv[i]);
    if (beta<thresh) continue;
    else{
      visindx[nvis]=i;
      vispid[nvis]=ipid;
      nvis++;
//      if (ipid<MAXNVIS) nvisarr[ipid]=nvisarr[ipid]+1;
      if ((ipid==5)||(ipid==6)) nvmu++;
      if ((ipid==2)||(ipid==3)) nve++;
      if ((ipid==9)||(ipid==8)) nvpip++;
      if ((ipid==11)||(ipid==12)) nvk++;
      if ((ipid==14)||(ipid==15)) nvp++;
      if ((ipid==7)) nvpi0++;
      if (ipid==1) nvgam++;
    }
  }
//  nvmu = nvisarr[5]+nvisarr[6];
//  nve  = nvisarr[2]+nvisarr[3];
//  nvpip = nvisarr[9]+nvisarr[8];
//  nvpi0 = nvisarr[7];
//  nvgam = nvisarr[1];
//  nvp   = nvisarr[14]+nvisarr[15];
//  nvk   = nvisarr[11]+nvisarr[12];
  nvoth=nvis-nvmu-nve-nvpip-nvpi0-nvgam-nvp-nvk;
  return;
}

float visRing::getbeta(int ipid, float pmom){
  if (ipid==0) return 0; //other thing?
//  if (ipid==4) return 0; //neurino 
//  if (ipid==14) return 0; //neutron 
  if (ipid==7) return 1; //pi0 always visible
  if ((ipid==1)&&(pmom>3.)) return 1;
  if ((ipid==1)&&(pmom<=3.)) return 0;
  float m = massof[ipid];
  //cout<<"mass: "<<m<<endl;
  float E = sqrt(m*m + pmom*pmom);
  //cout<<"beta: "<<pmom/E<<endl;
  return (pmom/E);
}

visRing::visRing(fqReader* fqin){
  FVflg = 0;
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

visRing::visRing(fqReaderFV* fqin){
  //set FV flg
  FVflg = 1;
  fqfv = fqin;
  massof[1] = 0;  //gamma
  massof[2] = 0.511e1; //positron
  massof[3] = 0.511e1; //electron 
  massof[4] = 0; //neutrino
  massof[5] = 0.1057e3; //mu+
  massof[6] = 0.1057e3; //mu-
  massof[7] = 0.1350e3; //pi0
  massof[8] = 0.1396e3; //pi+
  massof[9] = 0.1396e3; //pi-
  massof[10] = 0.4977e3; //K long
  massof[11] = 0.4937e3; //K+
  massof[12] = 0.4937e3; //K-
  massof[13] = 0.9396e3; //n
  massof[14] = 0.9383e3; //p
  massof[15] = 0.9383e3; //p-
  massof[17] = 0.5475e3; //eta
}

#endif
