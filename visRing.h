#ifndef VISRING_H
#define VISRING_H

#include "shared.h"
#include "fqEvent.h"
#ifdef T2K
#include "t2kfqEvent.h"
#endif
#include "TMath.h"
#include "TString.h"
#include "iostream"
#include "FVCalculators.cxx"
#include <map>

#define MAXNVIS 200

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
  int addvisible(int ipid, int index, double momentum,int flgscnd=0);
  void addvisiblesecondary(int ipid, int index, double momentum);
  void countdecaypi0();
  void initconstants();


 //useful quantities
 double getbeta(int ipid, double pmom);
 double getpcrit(int ipid);
 double getEnergy(int ipid, double pmom);
 double getEcrit(int ipid);
 double getKEcrit(int ipid);
 double getKE(int ipid, double pmom);
 double getVisibleEnergy(int ipid, double pmom);
 double getVisWall(int index, int scndflg);
 double getVisTowall(int index, int scndflg);
 void printsecondaryindex(int ipid);
 void printsecondaryinfo(int idx);
 int hasdschild(int vcindex);
 int pdg2geant(int ipid);
 map<int,double> massof;
 map<int,TString> nameof;
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
 double vistime[MAXNVIS]; //time of cretion of visible ring
 int visscndpid[MAXNVIS]; //pid of visible secondary ring
 int visscndparentid[MAXNVIS]; //parent pid of visible secondary ring
 int nvisscnd;  //number of visible secondary particles
 double visbrightness[MAXNVIS];  //strength of each visible ring
 double vismom[MAXNVIS]; //momentum of visible particle
 double viswall[MAXNVIS]; // true wall of particle that produced this ring
 double vistowall[MAXNVIS]; // true towall of particle that produced this ring
 double mumom[MAXNVIS]; //momentum of muons
 double pimom[MAXNVIS]; //momentum of charged pions
 double pi0mom[MAXNVIS]; //momentum of neutral pionsj
 double emom[MAXNVIS]; //momentum of e+-
 double gammom[MAXNVIS]; //momentum of gammas
 double protmom[MAXNVIS]; //momentum of proton
 double vismrbrightness; // multi-ring parameter (strength of 2nd strongest ring);
 double vismrt1; // time of MVR
 double vismrt2; // time of 2MVR
 int    vismrpid2; // pid of 2nd most visible ring
 int    vismrpid1; // pid of most visible ring
 int    vismrtype1; // 1 for shower, 0 for non-shower
 int    vismrtype2; // 1 for shower, 0 for non-shower
 double vismrwall1; 
 double vismrwall2;  
 double vismrtowall1;
 double vismrtowall2;
 double vismrwallmin;
 double vismrtowallmin;
 int nvisarr[MAXNVIS];

 // constants
 double Cthresh; //Cherenkov threshold in c
 double Tthresh; //cutoff time in ns to be counted in this event
 double gamthresh; // min energy of gamma
 int flgverbprime; //< verbose flag
 int flgverbscnd; //< verbose flag
 // debuggin
 void testevent(int iev);
};

#endif
