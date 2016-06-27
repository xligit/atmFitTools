#ifndef __COVXSEC_H__
#define __COVXSEC_H__

#include "covBase.h"

class covXsec : public covBase
{
 public:
  covXsec(std::string name, std::string file, unsigned int seed = 0);
  ~covXsec();

  double getLikelihood();
  void randomize();
  void correlateSteps();
  void proposeStep();
  void throwNominal(bool nomValues = true, unsigned int seed = 0);
  double GetWeightFrac(int i);
  double GetWeight(int i);

 protected:
  void InitPars(float stepScale, unsigned int seed);
  TVectorD *xsec_param_prior;
  TVectorD *xsec_param_nom;
  TVectorD *xsec_param_lb;
  TVectorD *xsec_param_ub;
  double *fParNom;
  double *fParType;

};

#endif
