#ifndef FVCALC_C
#define FVCALC_C

#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <iostream>

using namespace std;

//FV CALCULATORS//////////////////////////////////////////////
double calcWall(TVector3* vpos){
  double xx = vpos->X();
  double yy = vpos->Y();
  double zz = vpos->Z();
  double Rmax = 1690.;
  double Zmax = 1810.;
  double rr   = sqrt(xx*xx + yy*yy);
  double signflg = 1;
  if ((rr>Rmax)||(fabs(zz)>Zmax)) signflg  = -1;
  double wall = signflg*fmin(fabs(Zmax-fabs(zz)),fabs(Rmax-rr));
  if (wall<-3000.) return 3000.;
  return wall;
}



double calcToWallCustom(TVector3 *vpos, TVector3 *vdir, double dt){

  double R = 1690.0;
  double Z = 1810.0;
  double x0 = vpos->X();
  double y0 = vpos->Y();
  double z0 = vpos->Z();
  double dirx = vdir->X();
  double diry = vdir->Y();
  double dirz = vdir->Z();
  double r=0;
  double z=z0;
  double x=x0;
  double y=y0;
  double t=0.;
  int    flg =0;
  double tmax = 10000.;
  double sign = 1.;
  if (fabs(z0)>Z) return -10000;;
  if (sqrt(x0*x0+y0*y0)>R) return -10000;
  while ((flg==0)&&(t<tmax)){
    x = x0 + t*dirx;
    y = y0 + t*diry;
    z = z0 + t*dirz;
    r = sqrt(x*x + y*y);
    if (sign>0){
      if (fabs(z)>Z) flg=1;
      if (r>R) flg=1;
    }
    else{
      if (fabs(z)<Z) flg=1;
      if (r<R) flg=1;
    } 
    t+=dt;
  }
  //its the previous iteration that caused exit:
  t-=dt;
  x-=dt*dirx;
  y-=dt*diry;
  z-=dt*dirz;
  vpos->SetXYZ(x,y,z);

  return t*sign;
}

double calcPerimeter(TVector3 *vpos, TVector3* vdir){
  int npts = 20;
  TVector3 P[20]; //position vector for points
  for (int i=0;i<npts;i++){
     P[i].SetXYZ(vpos->X(),vpos->Y(),vpos->Z());
  }
  TVector3 D[20]; //direction vector for points
  TVector3 zhat;
  zhat.SetXYZ(0.,0.,1.);
  TVector3 perp; //vector perpendicular to direction
  perp = zhat.Cross(*vdir);
  double changle = 0.72245;
  double dangle = (2*3.14159)/(double)npts;
  double angle = 0.;
  TVector3 vdirtmp;
  vdirtmp.SetXYZ(vdir->X(),vdir->Y(),vdir->Z());
  TVector3 dir0 =vdirtmp;
  dir0.Rotate(changle,perp);
  for (int j=0;j<npts;j++){
    D[j] = dir0;
    D[j].Rotate(angle,vdirtmp);
    angle+=dangle;
  }
  double perim = 0;
  for (int k=0;k<npts;k++){
    calcToWallCustom(&P[k],&D[k],1.0); //sets P to value at ID wall
  }
  for (int pt=1;pt<npts;pt++){
    perim+=(P[pt-1]-P[pt]).Mag();
  }
  perim+=(P[npts-1]-P[npts]).Mag();
  return perim;
}

double calcToWall(TVector3* vpostmp, TVector3* vdirtmp){
  double towallval = 0.;
  TVector3 thepos;
  TVector3 thedir;
  thepos.SetXYZ(vpostmp->X(),vpostmp->Y(),vpostmp->Z());
  thedir.SetXYZ(vdirtmp->X(),vdirtmp->Y(),vdirtmp->Z());

  towallval+= calcToWallCustom(&thepos,&thedir,1.);
  towallval+= calcToWallCustom(&thepos,&thedir,0.1);
  return towallval;
}

double calcPhiWall(TVector3* vpos, TVector3* vdir){
  double Rmax = 1690.;
  double Zmax = 1810.;
  double sign = -1.;
  double R = sqrt((vpos->X()*vpos->X())+
                  (vpos->Y()*vpos->Y()));
  double Z = vpos->Z();
  TVector3 wallnorm;
  TVector3 rcdir;
  rcdir.SetXYZ(vdir->X(),vdir->Y(),vdir->Z());
  if (fabs(R-Rmax)>fabs(fabs(Z)-Zmax)) {
    //closer to endcaps
    if (Z<0.) wallnorm.SetXYZ(0.,0.,1);
    if (Z>=0.) wallnorm.SetXYZ(0.,0.,-1);
  }
  else {
    wallnorm.SetXYZ(-1*vpos->X(),-1*vpos->Y(),0.);
   // wallnorm*=sign;
  }
  return (180./TMath::Pi())*wallnorm.Angle(rcdir);  //return angle between norm and direction
}


//////////////////////////////////////////////////////////////


#endif
