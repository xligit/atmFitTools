#ifndef FVCALC_C
#define FVCALC_C

#include "TPolyLine3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <iostream>
#include "TRandom2.h"
#include "TCanvas.h"

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

double calcWall2(TVector3* vpos){
  double xx = vpos->X();
  double yy = vpos->Y();
  double zz = vpos->Z();
  double Rmax = 1690.;
  double Zmax = 1810.;
  double rr   = sqrt(xx*xx + yy*yy);
  double absz = TMath::Abs(zz);
  //check if vertex is outside tank
  double signflg = 1.;
  if (absz>Zmax) signflg = -1.;
  if (rr>Rmax)   signflg = -1.;
  //find min distance to wall
  double distz = TMath::Abs(Zmax-absz);
  double distr = TMath::Abs(Rmax-rr);
  double wall = signflg*fmin(distz,distr);
  return wall;
}

double calcToWallCustom(TVector3 *vpos, TVector3 *vdir, double dt){

  double R = 1690.0;
  double Z = 1810.0;
  double x0 = vpos->X();
  double y0 = vpos->Y();
  double z0 = vpos->Z();
  double r0 = sqrt( x0*x0 + y0*y0);
  double dirx = vdir->X();
  double diry = vdir->Y();
  double dirz = vdir->Z();
  double r=0;
  double z=0;
  double x=0;
  double y=0;
  double t=0.;
  int    flg =0;
  double tmax = 100000.;
  if (fabs(z0)>Z) return -10000;
  if (r0>R) return -10000;
  while ((flg==0)&&(t<tmax)){
    x = x0 + t*dirx;
    y = y0 + t*diry;
    z = z0 + t*dirz;
    r = sqrt( (x*x) + (y*y) );
    if (fabs(z)>Z) flg=1;
    if (r>R) flg=1;
    t+=dt;
  }
  //its the previous iteration that caused exit:
  t-=dt;
  x-=dt*dirx;
  y-=dt*diry;
  z-=dt*dirz;
  vpos->SetXYZ(x,y,z);
  return t;
}





/////////////////////////////////////////////////////
// estimate the perimeter of the cherenkov ring
double calcMinCone(TVector3 *vpos, TVector3* vdir, int npts=50){

  //////////////////////////////////////////////////////////////
  // each ray starts at initial point defined from given vertex
  TVector3 P[100]; // initial position vector for rays
  for (int i=0;i<npts;i++){
     P[i].SetXYZ(vpos->X(),vpos->Y(),vpos->Z());
  }


  ////////////////////////////////////////////////////////////
  // set the ray directions
  // we must first find a vector that is perpendicular to particle direction
  int bestaxis=0;
  if (TMath::Abs(vdir->Y())<TMath::Abs(vdir->X())){
    bestaxis=1;
    if (TMath::Abs(vdir->Z())<TMath::Abs(vdir->Y())) bestaxis=2;
  }
  else{
    bestaxis = 0;
    if (TMath::Abs(vdir->Z())<TMath::Abs(vdir->X())) bestaxis=2;
  }
  TVector3 vhat;
  if (bestaxis==0)vhat.SetXYZ(1.,0.,0.);
  if (bestaxis==1)vhat.SetXYZ(0.,1.,0.);
  if (bestaxis==2)vhat.SetXYZ(0.,0.,1.);
  TVector3 vperp = vhat.Cross(*vdir); // this vector is perpendicular to particle direction

  // now we make initial ray, and then rotate it around particle direction
  TVector3 D[100]; //direction vector for points
  double changle = 0.72245; // Cherenkov angle in water
  double dangle = (2*3.14159)/(double)npts; // angular differnce between rays
  double angle = 0.;  // inital angle
  TVector3 dir0; // inital ray
  dir0.SetXYZ(vdir->X(),vdir->Y(),vdir->Z()); // start equalt to particle dir
  dir0.Rotate(changle,vperp); // rotate it onto Cherenkov cone
  for (int j=0;j<npts;j++){
    D[j] = dir0;  // put direction on cone
    D[j].Rotate(angle,*vdir); // rotate on cone
    angle+=dangle;
  }


  //////////////////////////////////////////////////
  // now that rays are defined, let's calculate the perimeter
  double raymin = 10000;
  for (int k=0;k<npts;k++){
    // calc dist to wall within 1 cm;	  
    double raytowall = calcToWallCustom(&P[k],&D[k],1.0); //sets P to value at ID wall
    if (raytowall<raymin) raymin=raytowall;
  }


  return raymin;
}





/////////////////////////////////////////////////////
// estimate the perimeter of the cherenkov ring
double calcPerimeter(TVector3 *vpos, TVector3* vdir, int npts=50, int visflg=0){
  

  //////////////////////////////////////////////////////////////
  // each ray starts at initial point defined from given vertex
  TVector3 P[100]; // initial position vector for rays
  for (int i=0;i<npts;i++){
     P[i].SetXYZ(vpos->X(),vpos->Y(),vpos->Z());
  }

  ////////////////////////////////////////////////////////////
  // set the ray directions

  // we must first find a vector that is perpendicular to particle direction
  int bestaxis=0;
  if (TMath::Abs(vdir->Y())<TMath::Abs(vdir->X())){
    bestaxis=1;
    if (TMath::Abs(vdir->Z())<TMath::Abs(vdir->Y())) bestaxis=2;
  }
  else{
    bestaxis = 0;
    if (TMath::Abs(vdir->Z())<TMath::Abs(vdir->X())) bestaxis=2;
  }
  TVector3 vhat;
  if (bestaxis==0)vhat.SetXYZ(1.,0.,0.);
  if (bestaxis==1)vhat.SetXYZ(0.,1.,0.);
  if (bestaxis==2)vhat.SetXYZ(0.,0.,1.);
  TVector3 vperp = vhat.Cross(*vdir); // this vector is perpendicular to particle direction

  // now we make initial ray, and then rotate it around particle direction
  TVector3 D[100]; //direction vector for points
  double changle = 0.72245; // Cherenkov angle in water
  double dangle = (2*3.14159)/(double)npts; // angular differnce between rays
  double angle = 0.;  // inital angle
  TVector3 dir0; // inital ray
  dir0.SetXYZ(vdir->X(),vdir->Y(),vdir->Z()); // start equalt to particle dir
  dir0.Rotate(changle,vperp); // rotate it onto Cherenkov cone
  for (int j=0;j<npts;j++){
    D[j] = dir0;  // put direction on cone
    D[j].Rotate(angle,*vdir); // rotate on cone
    angle+=dangle;
  }

  //////////////////////////////////////////////////
  // now that rays are defined, let's calculate the perimeter
  double perim = 0;
  for (int k=0;k<npts;k++){
    // calc dist to wall within 1 cm;	  
    calcToWallCustom(&P[k],&D[k],1.0); //sets P to value at ID wall
  }
  //loop over points and add distances
  for (int pt=1;pt<npts;pt++){
    perim+=(P[pt-1]-P[pt]).Mag();
  }
  //don't forget to connect first  point to last
  perim+=(P[0]-P[npts-1]).Mag();


  // draw it if visflg (for debugging)
  if (visflg){
    
   //make arrays for polylines
   double x[2];
   double y[2];
   double z[2];
   TPolyLine3D *rays[100];
   for (int i=0; i<npts; i++){
     x[0] = vpos->X();
     y[0] = vpos->Y();
     z[0] = vpos->Z();
     x[1] = P[i].X();
     y[1] = P[i].Y();
     z[1] = P[i].Z();
     rays[i] = new TPolyLine3D(2,x,y,z);
     if (i==0) rays[0]->Draw();
     else{
       rays[i]->Draw("same");
     }
   }
  }

  // return final perimeter value

  return perim;
}

double calcToWall(TVector3* vpostmp, TVector3* vdirtmp){
  double towallval = 0.;
  TVector3 thepos;
  TVector3 thedir;
  thepos.SetXYZ(vpostmp->X(),vpostmp->Y(),vpostmp->Z());
  thedir.SetXYZ(vdirtmp->X(),vdirtmp->Y(),vdirtmp->Z());
  double wallval = calcWall(vpostmp)
  if (wallval<0) return wallval; 
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


void towallValidataion(int nrays){

  TCanvas* cc = new TCanvas("cc","cc",800,700);

  const int NRAYS = nrays;

  TPolyLine3D *rays[NRAYS];

  TVector3 vpos(0,0,0);
  TVector3 *vdir[NRAYS];

   TRandom2 * randy = new TRandom2(nrays);

  //make random directions
  for (int i=0; i<nrays; i++){
    double xdir = randy->Gaus(0,1);
    double ydir = randy->Gaus(0,1);
    double zdir = randy->Gaus(0,1);
    vdir[i] = new TVector3(xdir,ydir,zdir);
    vdir[i]->SetMag(1.);

    //calc towall
    TVector3* vposf = new TVector3(vpos.X(),vpos.Y(),vpos.Z());
    calcToWallCustom(vposf,vdir[i],1.);
    //make line
    double xx[2];
    double yy[2];
    double zz[2];
    xx[0] = vpos.X();
    yy[0] = vpos.Y();
    zz[0] = vpos.Z();
    xx[1] = vposf->X();
    yy[1] = vposf->Y();
    zz[1] = vposf->Z();
    rays[i] = new TPolyLine3D(2,xx,yy,zz);
    if (i==0) rays[0]->Draw();
    else{
      rays[i]->Draw("same");
    }
  }

  cc->Print("~/fvtest.png");

  return;

}

#endif
