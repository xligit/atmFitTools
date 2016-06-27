#ifndef _COVBASE_H_
#define _COVBASE_H_

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>

#include "TFile.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"

class covBase 
{
 public:
  covBase();
  covBase(std::string name, std::string file);
  covBase(std::string name, std::string file, unsigned int seed);
  virtual ~covBase();

  void setName(std::string name) { mName = name; }
  void setParName(int i, std::string name) { fParName[i] = name; }
  void setEvalLikelihood(int i, bool e) { fParEvalLikelihood[i] = e; }
  void setPar(int i, double val) { fParCurr[i] = val; }
  void setPars(std::vector<double> pars);
  void setThrowPar(int i, double val) { fParSigma[i] = val; } // set sigma of a parameter or decide whether to throw a parameter or not
  void setStepScales(float scale);
  void setStepScale(int i, float scale) { fStepScale[i] = scale; }
  void setPropFunc(int i, TF1 *func);

  TMatrixDSym* getCovMatrix() { return cov; }
  TMatrixDSym* getInvCovMatrix() { return invCov; }
  std::string getName() { return mName; }
  std::string getParName(int i) { return fParName[i]; }
  std::vector<double> getNominalArray();
  std::vector<double> getProposedArray();
  std::vector<double> getCurrentArray();
  double getNominal(int i) { return nominal.at(i); }
  double getUncertainty(int i) {return fParSigma[i]; }
  double getUp(int i) { return fParUp[i]; }
  double getLow(int i) { return fParLow[i]; }
  double getProposed(int i) { return fParProp[i]; }
  double getCurrent(int i) { return fParCurr[i]; }
  double getInit(int i) {return fParInit[i]; }
  int getNPar() { return size; }
  int getNXsecPar() {return size;}
  TF1* getPropFunc(int i) { return fPropKernel[i]; }
  bool isParameterFixed(int i) { return fParEvalLikelihood[i]; }

  void PrintNominal();
  void PrintPars();

  virtual void throwNominal(bool nomValues = true, unsigned int seed = 0);
  virtual double getLikelihood();
  virtual void proposeStep();
  virtual double GetWeightFrac(int i) {return fParCurr[i]-1;}
  virtual double GetWeight(int i) {return fParCurr[i];}
  void acceptStep();

 protected:
  void init(std::string name, std::string file);
  void correlateSteps();
  void genPropKernels();
  void CholeskyDecomp(int npars, TMatrixD &chel_matrix);
  void randomize();

  TRandom3 *rnd;
  int size;
  int size_xsec;
  std::string mName;
  TMatrixDSym *cov;
  TMatrixDSym *invCov;
  std::vector<double> nominal;

  TVectorD randPar;
  TMatrixD *chel;
  float *fStepScale;
  TF1 **fPropKernel;
  std::string *fParName;
  double *fParInit;
  double *fParCurr;
  double *fParProp;
  double *fParSigma;
  double *fParUp;
  double *fParLow;
  bool *fParEvalLikelihood;

};

#endif
