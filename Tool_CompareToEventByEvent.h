#include "histoManager.h"
#include "hSplines.h"
#include "atmFitPars.h"

using namespace std;

class compareToEventByEvent{
  public:
  compareToEventByEvent(atmFitPars* fitpars, TChain* tr, const char* hfilename, const char* sfilename, bool separateNeutMode = false);
  histoManager* hManager;
  TH1D* htmp1;
  TH1D* htmp2;
  TH1D* htmp3;
  TChain* mcTree;
  atmFitPars* thePars;
  bool separateNeutMode;
  fQreader* mcEvt;
  float getEvtWeight(int isys);
  float att[1000];
  void fillAttributes();
  void comparePrediction(int isamp,int ibin, int icomp, int iatt, int isys);
  void comparePrediction(int isamp,int ibin, int icomp, int imode, int iatt, int isys); 
  void comparePrediction(int isamp,int ibin, int icomp, int imode, int iatt, bool all);
};

