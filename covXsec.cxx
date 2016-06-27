#include "TObjArray.h"
#include "TObjString.h"
#include "TDecompChol.h"

#include "ThrowParms.h"
#include "covXsec.h"

covXsec::covXsec(std::string name, std::string file, unsigned seed)
  : covBase(name, file, seed)
{
  TFile infile(file.c_str(), "read");
  xsec_param_nom = (TVectorD*)infile.Get("xsec_param_nom");
  xsec_param_prior = (TVectorD*)infile.Get("xsec_param_prior");
  //xsec_param_id = (TMatrixD*)infile.Get("xsec_param_id");
  xsec_param_lb = (TVectorD*)infile.Get("xsec_param_lb");
  xsec_param_ub = (TVectorD*)infile.Get("xsec_param_ub");
  TObjArray * arr = (TObjArray*)infile.Get("xsec_param_names");
  fParType = new double[size];
  fParNom = new double[size];
  for (int i = 0; i < size; ++i) {
    fParName[i] = std::string(((TObjString*)arr->At(i))->GetName());
    fParUp[i] = (*xsec_param_ub)(i);
    fParLow[i] = (*xsec_param_lb)(i);
    nominal.at(i) = (*xsec_param_prior)(i);
    fParNom[i] = (*xsec_param_nom)(i);
    if (((TObjString*)arr->At(i))->GetString().Contains("RPA")) fParType[i] = 2;
    else fParType[i] = 0;
    if (fParNom[i] != 0) {
      (*xsec_param_prior)(i) /= fParNom[i];
      (*xsec_param_lb)(i) /= fParNom[i];
      (*xsec_param_ub)(i) /= fParNom[i];
      fParUp[i] /= fParNom[i];
      fParLow[i] /= fParNom[i];
      fParNom[i] = 1;
    }
  }
  InitPars(0.001, seed);
  std::cout<<"----------- nominal (unscaled) ------------\n";
  for (int i = 0; i < size; ++i) {
    std::cout << i << " " << fParName[i] << " = " << (*xsec_param_nom)(i) << ", " << sqrt((*cov)(i,i)) << "\n";
  }
  std::cout<<"--------------------------------"<<std::endl;

}

covXsec::~covXsec()
{
}

void covXsec::InitPars(float stepScale, unsigned int seed)
{
  for (int i = 0; i < size; ++i) {
    fParInit[i] = (*xsec_param_prior)(i);
    fParCurr[i] = fParInit[i];
    fParProp[i] = fParCurr[i];
    //fParSigma[i] = stepScale;
    fParSigma[i] = TMath::Sqrt((*cov)(i,i));
    fStepScale[i] = stepScale;
  }
  //genPropKernels();
  throwNominal(true, seed);
  randomize();
  correlateSteps();
}

double covXsec::getLikelihood()
{
  double LogL = 0;
  for (int i = 0; i < size; ++i) {
    if (fParProp[i] > fParUp[i] || fParProp[i] < fParLow[i]) return 99999.9;
    for (int j = 0; j < size; ++j) {
      if (fParEvalLikelihood[i] && fParEvalLikelihood[j] && fParType[i]==0 && fParType[j]==0) {
	LogL += (nominal.at(i)-fParProp[i]) * (nominal.at(i)-fParProp[i]) * (*invCov)(i,j) * 0.5;
      }
    }
  }
  return LogL;
}

void covXsec::proposeStep()
{
  randomize();
  correlateSteps();
}

void covXsec::randomize()
{
  for (int i = 0; i < size; ++i) {
    if (fParSigma[i] > 0 && fParType[i]==0) randPar(i) = rnd->Gaus(0, 1); // fPropKernel[i]->GetRandom();
    else randPar(i) = 0;
  }
}

void covXsec::correlateSteps()
{
  TVectorD corr_throw = (*chel) * randPar;
  for (int i = 0; i < size; ++i) {
    if (fParSigma[i] > 0 && fParType[i]==0) {
      fParProp[i] = fParCurr[i] + corr_throw(i) * fStepScale[i];
    } else if (fParSigma[i] > 0 && fParType[i]==2) {
      double rn = rnd->Uniform(0, 1);
      if (rn < 0.5) {
	if (fParCurr[i]==0) fParProp[i] = 1;
	else if (fParCurr[i]==1) fParProp[i] = 0;
      } else {
	fParProp[i] = fParCurr[i];
      }
    } else {
      fParProp[i] = fParCurr[i];
    } 
  }
}

void covXsec::throwNominal(bool nomValues, unsigned int seed)
{
  TDecompChol chdcmp(*cov);
  if (!chdcmp.Decompose()) {
    std::cerr << "Cholesky decomposition failed for "<< mName << std::endl;
    exit(-1);
  }
  chel = new TMatrixD(chdcmp.GetU());
  CholeskyDecomp(size, *chel);
  ThrowParms *nom_throws = new ThrowParms(*xsec_param_prior, *cov);
  nom_throws->SetSeed(seed);
  nominal.clear();
  nominal.resize(size);
  if (!nomValues) {
    bool throw_again = true;
    while (throw_again) {
      throw_again = false;
      nom_throws->ThrowSet(nominal);
      for (int i = 0; i < size; ++i) {
	if (fParSigma[i] < 0) {
	  nominal.at(i) = (*xsec_param_prior)(i);
	  continue;
	}
	if (nominal.at(i) < fParLow[i]) throw_again = true;
	if (nominal.at(i) > fParUp[i])  throw_again = true;
      }
    }
  } else {
    for (int i = 0; i < size; ++i) nominal.at(i) = (*xsec_param_prior)(i);
  }
}

double covXsec::GetWeightFrac(int i)
{
  if (i >= 0 && i <= size) return fParProp[i] - 1.;
  else return 1;
}

double covXsec::GetWeight(int i)
{
  if (i >= 0 && i <= size) return fParProp[i];
  else return 1;
}
