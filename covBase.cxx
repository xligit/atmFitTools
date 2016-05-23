#include <math.h>
#include "TDecompChol.h"
#include "ThrowParms.h"
#include "covBase.h"

covBase::covBase() {}

covBase::covBase(std::string name, std::string file)
{
  rnd = new TRandom3(1);
  init(name, file);
}

covBase::covBase(std::string name, std::string file, unsigned int seed)
{
  rnd = new TRandom3(seed);
  init(name, file);
}

covBase::~covBase()
{
  delete fParName;
  delete fParInit;
  delete fParSigma;
  delete fParCurr;
  delete fParProp;
  delete fParUp;
  delete fParLow;
  delete fParEvalLikelihood;
  delete fStepScale;
  delete[] fPropKernel;

}

void covBase::init(std::string name, std::string file)
{
  TFile infile(file.c_str(), "read");
  TMatrixDSym *covMatrix = (TMatrixDSym*)(infile.Get(name.c_str())); 
  cov = (TMatrixDSym*)covMatrix->Clone(name.c_str());
  invCov = (TMatrixDSym*)covMatrix->Clone("invCov");
  invCov->Invert();
  mName = name;
  size = cov->GetNrows();
  size_xsec = size;
  fParName = new std::string[size];
  fParInit = new double[size];
  fParSigma = new double[size];
  fParCurr = new double[size];
  fParProp = new double[size];
  fParUp = new double[size];
  fParLow = new double[size];
  fParEvalLikelihood = new bool[size];
  fPropKernel = new TF1*[size];
  fStepScale = new float[size];

  for (int i = 0; i < size; ++i) {
    nominal.push_back(1.);
    fParEvalLikelihood[i] = true;
    fStepScale[i] = 1.;
  }

  randPar.ResizeTo(size);

}

void covBase::throwNominal(bool nomValues, unsigned int seed)
{
  TDecompChol chdcmp(*cov);
  if (!chdcmp.Decompose()) {
    std::cerr << "Cholesky decomposition failed for " << mName << std::endl;
    exit(-1);
  }

  chel = new TMatrixD(chdcmp.GetU());
  CholeskyDecomp(size, *chel);
  TVectorD vec(size);
  for (int i = 0; i < size; ++i) vec(i) = 1.;

  ThrowParms *nom_throws = new ThrowParms(vec, (*cov));
  nom_throws->SetSeed(seed);
  nominal.clear();
  nominal.resize(size);

  if (!nomValues) {
    bool throw_again = true;
    while (throw_again) {
      throw_again = false;
      nom_throws->ThrowSet(nominal);
      for (int i = 0; i < size; ++i) {
	if (fParSigma[i] < 0) { // if parameter is fixed, do not throw
	  nominal.at(i) = 1.0;
	  continue;
	} 
	if (nominal.at(i) < 0) {
	  nominal.at(i) = 0;
	  throw_again = true;
	}
      }
    }
  }
  delete nom_throws;
}

void covBase::CholeskyDecomp(int npars, TMatrixD &chel_mat)
{
  for(int i=0; i<npars; i++)
    for(int j=0; j<npars; j++){
      //if diagonal element  
      if(i==j){
        chel_mat(i,i) = (*cov)(i,i);
        for(int k=0; k<=i-1; k++) chel_mat(i,i) = chel_mat(i,i)-pow(chel_mat(i,k),2);
        chel_mat(i,i) = sqrt(chel_mat(i,i));
        //if lower half 
      } else if(j<i) {
        chel_mat(i,j) = (*cov)(i,j);
        for(int k=0; k<=j-1; k++) chel_mat(i,j) = chel_mat(i,j)-chel_mat(i,k)*chel_mat(j,k);
        chel_mat(i,j) = chel_mat(i,j)/chel_mat(j,j);
      } else chel_mat(i,j) = 0.;
    }
}

void covBase::genPropKernels()
{
  for (int i = 0; i < size; ++i) {
    fPropKernel[i] = NULL;
    if (fParSigma[i] > 0) {
      double e = sqrt((*cov)(i,i));
      fPropKernel[i] = new TF1(Form("%s_%03d", mName.c_str(), i), "gaus", -e*10., e*10.);
      fPropKernel[i]->SetParameters(1./(e*sqrt(2*3.14159265)), 0., e);
    }
  }
}

void covBase::setPropFunc(int i, TF1 *func)
{
  if (i >= size) {
    std::cerr << "Parameter index out of range." << std::endl;
    exit(-1);
  }
  fPropKernel[i] = func;
  fParSigma[i] = func->GetParameter(1);
}

void covBase::proposeStep()
{
  randomize();
  correlateSteps();
}

void covBase::randomize()
{
  for (int i = 0; i < size; ++i) {
    if (fParSigma[i] > 0) {
      //randPar(i) = rnd->Gaus(0,1);
      randPar(i) = fPropKernel[i]->GetRandom();
    } else {
      randPar(i) = 0;
    }
  }
} 

double covBase::getLikelihood()
{
  double LogL = 0;
  for (int i = 0; i < size; ++i) {
    if (fParProp[i] < 0) return 99999.9;
    for (int j = 0; j < size; ++j) {
      if (fParEvalLikelihood[i] && fParEvalLikelihood[j]) LogL += (fParProp[i] - nominal[i]) * (fParProp[j] - nominal[j]) * (*invCov)(i,j) * 0.5;
    }
  }
  return LogL;
}

void covBase::correlateSteps()
{
  TVectorD corr_throw = (*chel)*randPar;
  for (int i = 0; i < size; ++i) {
    if (fParSigma[i] > 0) fParProp[i] = fParCurr[i] + corr_throw(i)*fStepScale[i];
  }
}

void covBase::acceptStep()
{
  for (int i = 0; i < size; ++i) fParCurr[i] = fParProp[i];
}

std::vector<double> covBase::getNominalArray()
{
  return nominal;
}

std::vector<double> covBase::getProposedArray()
{
  std::vector<double> prop;
  for (int i = 0; i < size; ++i) prop.push_back(fParProp[i]);
  return prop;
}

std::vector<double> covBase::getCurrentArray()
{
  std::vector<double> curr;
  for (int i = 0; i < size; ++i) curr.push_back(fParCurr[i]);
  return curr;
}

void covBase::setPars(std::vector<double> pars)
{
  if (pars.size() != (unsigned int)size) {
    std::cerr << "Error with covBase::setPars(std::vector<double> pars): pars.size() != size."<<std::endl;
    exit(-1);
  }
  for (int i = 0; i < size; ++i) fParCurr[i] = pars.at(i);
}

void covBase::setStepScales(float scale)
{
  for (int i = 0; i < size; ++i) fStepScale[i] = scale;
}

void covBase::PrintNominal()
{
  std::cout<<"----------- Nominal ------------\n";
  for (int i = 0; i < size; ++i) {
    std::cout << i << " " << fParName[i] << " = " << nominal.at(i) << ", sigma = " << fParSigma[i] << ", low = " << fParLow[i] << ", up = " << fParUp[i] << "\n";
  }
  std::cout<<"--------------------------------"<<std::endl;
}

void covBase::PrintPars()
{
  std::cout<<"----------- Current ------------\n";
  for (int i = 0; i < size; ++i) {
    std::cout << i << " " << fParName[i] << " = " << fParCurr[i] << "\n";
  }
  std::cout<<"--------------------------------"<<std::endl;
}
