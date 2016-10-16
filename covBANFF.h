#ifndef __COVBANFF_H__
#define __COVBANFF_H__

#include "TAxis.h"
#include "covBase.h"

class covBANFF : public covBase
{
 public:
  covBANFF(std::string name, std::string file, unsigned int seed = 0, bool postfit = true);
  ~covBANFF();

  double getLikelihood();
  void randomize();
  void correlateSteps();
  void proposeStep();
  void throwNominal(bool nomValues = true, unsigned int seed = 0);
  double GetWeightFrac(int i);
  double GetWeight(int i);
  double getPrior(int i);

 protected:
  void InitPars(float stepScale, unsigned int seed);

  TVectorD *banff_param_pre;
  TVectorD *banff_param_prior;
  TVectorD *banff_param_nom;
  TVectorD *banff_param_lb;
  TVectorD *banff_param_ub;
  TAxis *nu_numubins;
  TAxis *nu_numubarbins;
  TAxis *nu_nuebins;
  TAxis *nu_nuebarbins;
  TAxis *antinu_numubins;
  TAxis *antinu_numubarbins;
  TAxis *antinu_nuebins;
  TAxis *antinu_nuebarbins;
  TAxis *cccohbins;
  TAxis *nccohbins;
  TAxis *ncotherbins;
  TAxis *nc1gammabins;
  TString *banffnames;

  double linex[2];
  double *fParPrefit;
  int *fParType;
  double *fParPrior;
  int nBeamPars;
};

#endif
