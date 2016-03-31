#ifndef __VISRING_H__
#define __VISRING_H__

#include "fqReader.h"
#include "fqReaderFV.h"
#include "TMath.h"
#include "iostream"
#include <map>

#define MAXNVIS 100

using namespace std;

//a class for counting the number of visible rings in an event

class visRing{
 public:
 //constructor//
 //construct a visring object using a fitqun reader object
 visRing(fqReader* fqin);
 visRing(fqReaderFV* fqin);
 //
 fqReader* fq;
 fqReaderFV* fqfv;
 void fillVisVar(); //fills the variables relating to number of visible
 //rings
 void fillVisVarFV(); //fills the variables relating to number of visible


 //flag for using FV tree
 int FVflg;

 //useful quantities
 float getbeta(int ipid, float pmom);
 map<int,float> massof;

 //visible variables
 int nvis;  //total number of visible rings
 int nvmu;  //number of visible muon rings
 int nve;   //number of visible electron rings
 int nvpip; //number of visible charged pion rings
 int nvpi0; //number of visible neutral pion rings
 int nvoth; //number of other visible rings
 int nvgam; //number of visible gamma rings
 int nvp;   //number of visible proton rings
 int nvk;   //numbr of visible kaon rings
 int nvlam; //number of lambda particle
 int visindx[MAXNVIS]; //index of each visible ring in particle stack
 int vispid[MAXNVIS]; //pid code of each visible ring
 int nvisarr[MAXNVIS];
 
};

#endif
