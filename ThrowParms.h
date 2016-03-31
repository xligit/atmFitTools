#ifndef __THROWPARAM_H__
#define __THROWPARAM_H__

#include <iostream>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <vector>

#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TDecompChol.h"

#include "TRandom3.h"
#include "TF1.h"

class ThrowParms 
{
 public:
  ThrowParms(TVectorD &parms, TMatrixDSym &covm);
  ~ThrowParms();   
  void SetSeed(int seed = 0) {rand.SetSeed(seed);};
  int GetSize() {return npars;};
  void ThrowSet(std::vector<double> &parms);

 private:
  void CheloskyDecomp(TMatrixD &chel_mat);
  void StdNormRand(double *z);
  //typedef TVectorT<double> TVectorD;
  //typedef TMatrixTSym<double> TMatrixDSym;
  //typedef TMatrixT<double> TMatrixD;
  TVectorD    *pvals;
  TMatrixDSym *covar; 
  TMatrixD    *chel_dec;
  int         npars;
  TRandom3    rand; 
  TF1         *gauss;
};
#endif
