#ifndef VISRING_H
#define VISRING_H

#include "shared.h"
#include "fqEvent.h"
#ifdef T2K
#include "t2kfqEvent.h"
#endif
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
#ifdef T2K
  visRing(t2kfqEvent *fqin);
#else
  visRing(fqEvent* fqin);
#endif

#ifdef T2K
  t2kfqEvent *fq;
#else
  fqEvent* fq;
#endif
//  visRing(TTree* tr);
  void fillVisVar(); //fills the variables relating to number of visible
  void countsecondary();
  void countprimary();
  void countprimaryvc();
  void calcderived();
  void addvisible(int ipid, int index, double momentum);
  void addvisiblesecondary(int ipid, int index, double momentum,int flgverb=0);
  void countdecaypi0();
  void initconstants();


 //useful quantities
 double getbeta(int ipid, double pmom);
 double getpcrit(int ipid);
 int pdg2geant(int ipid);
 map<int,double> massof;
 double showerthresh;

 //visible variables
 int nvis;  //total number of visible rings
 int nvisprime; // visible rings in primary stack
 int nvmu;  //number of visible muon rings
 int nve;   //number of visible electron rings
 int nvpip; //number of visible charged pion rings
 int nvpi0; //number of visible neutral pion rings
 int nvoth; //number of other visible rings
 int nvgam; //number of visible gamma rings (EM showers)
 int nvp;   //number of visible proton rings
 int nvk;   //numbr of visible kaon rings
 int visindx[MAXNVIS]; //index of each visible ring in particle stack
 int vispid[MAXNVIS]; //pid code of each visible ring
 int visscndpid[MAXNVIS]; //pid of visible secondary ring
 int visscndparentid[MAXNVIS]; //parent pid of visible secondary ring
 int nvisscnd;  //number of visible secondary particles
 double visstr[MAXNVIS];  //strength of each visible ring
 double vismom[MAXNVIS]; //momentum of visible particle
 double mumom[MAXNVIS];
 double pimom[MAXNVIS];
 double pi0mom[MAXNVIS];
 double emom[MAXNVIS];
 double gammom[MAXNVIS];
 double protmom[MAXNVIS];
 double vismrpar; // multi-ring parameter (strength of 2nd strongest ring);
 int    vismrpid2; // pid of 2nd most visible ring
 int    vismrpid1; // pid of most visible ring
 int nvisarr[MAXNVIS];

 // constants
 double Cthresh; //Cherenkov threshold in c
 double Tthresh; //cutoff time in ns to be counted in this event
 double gamthresh; // min energy of gamma

 // debuggin
 void testevent(int iev);
};

#endif
