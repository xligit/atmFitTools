#include "TObjArray.h"
#include "TObjString.h"
#include "TDecompChol.h"

#include "ThrowParms.h"
#include "covBANFF.h"

covBANFF::covBANFF(std::string name, std::string file, unsigned seed, bool postfit)
  : covBase(name, file, seed)
{
  TFile infile(file.c_str(), "read"); 

  if (postfit) {
    banff_param_nom = (TVectorD*)infile.Get("postfit_params")->Clone();
    banff_param_prior = (TVectorD*)infile.Get("postfit_params")->Clone();
    banff_param_pre = (TVectorD*)infile.Get("prefit_params")->Clone();
    std::cout<<"Using post-fit parameters"<<std::endl;
  } else {
    banff_param_nom = (TVectorD*)infile.Get("prefit_params")->Clone();
    banff_param_prior = (TVectorD*)infile.Get("prefit_params")->Clone();
    banff_param_pre = (TVectorD*)infile.Get("prefit_params")->Clone();
    std::cout<<"Using pre-fit parameters"<<std::endl;
  }
  
  nu_numubins = (TAxis*)infile.Get("sk_numode_numu_bins");
  nu_numubarbins = (TAxis*)infile.Get("sk_numode_numub_bins");
  nu_nuebins = (TAxis*)infile.Get("sk_numode_nue_bins");
  nu_nuebarbins = (TAxis*)infile.Get("sk_numode_nueb_bins");
  antinu_numubins = (TAxis*)infile.Get("sk_anumode_numu_bins");
  antinu_numubarbins = (TAxis*)infile.Get("sk_anumode_numub_bins");
  antinu_nuebins = (TAxis*)infile.Get("sk_anumode_nue_bins");
  antinu_nuebarbins = (TAxis*)infile.Get("sk_anumode_nueb_bins");
  cccohbins = (TAxis*)infile.Get("CCCOH_O_bins");
  nccohbins = (TAxis*)infile.Get("NCCOH_bins");
  ncotherbins = (TAxis*)infile.Get("NCOTHER_bins");
  nc1gammabins = (TAxis*)infile.Get("NC1GAMMA_bins");

  const int nbanffpars = size;
  banffnames = new TString[nbanffpars];
  int index = 0;
  // Lower bounds and upper bounds have to be hard coded because they are not in the BANFF file
  // It sucks for no good reason
  for (int i = 0; i < nu_numubins->GetNbins(); ++i) {
    banffnames[i+index] = Form("FHC #nu_{#mu} %d", i);
    setParName(i+index, Form("fhc_numu_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += nu_numubins->GetNbins();
  for (int i = 0; i < nu_numubarbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("FHC #bar{#nu}_{#mu} %d", i);
    setParName(i+index, Form("fhc_numub_%d", i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += nu_numubarbins->GetNbins();
  for (int i = 0; i < nu_nuebins->GetNbins(); ++i) {
    banffnames[i+index] = Form("FHC #nu_{e} %d", i);
    setParName(i+index, Form("fhc_nue_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += nu_nuebins->GetNbins();
  for (int i = 0; i < nu_nuebarbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("FHC #bar{#nu}_{e} %d", i+1);
    setParName(i+index, Form("fhc_nueb_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += nu_nuebarbins->GetNbins();
  linex[0] = index; //------------------------------------------------                       
  std::cout<<"there are "<<linex[0]<<" nu-mode flux parameters."<<std::endl;
  for (int i = 0; i < antinu_numubins->GetNbins(); ++i) {
    banffnames[i+index] = Form("RHC #nu_{#mu} %d", i+1);
    setParName(i+index, Form("rhc_numu_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += antinu_numubins->GetNbins();
  for (int i = 0; i < antinu_numubarbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("RHC #bar{#nu}_{#mu} %d", i+1);
    setParName(i+index, Form("rhc_numu_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += antinu_numubarbins->GetNbins();
  for (int i = 0; i < antinu_nuebins->GetNbins(); ++i) {
    banffnames[i+index] = Form("RHC #nu_{e} %d", i+1);
    setParName(i+index, Form("rhc_nue_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += antinu_nuebins->GetNbins();
  for (int i = 0; i < antinu_nuebarbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("RHC #bar{#nu}_{e} %d", i+1);
    setParName(i+index, Form("rhc_nueb_%d",i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += antinu_nuebarbins->GetNbins();
  nBeamPars = index + 1;
  linex[1] = index; //------------------------------------------------               
  std::cout<<"there are "<<(linex[1]-linex[0])<<" antinu-mode flux parameters."<<std::endl;
  // MAQE
  banffnames[index] = "MAQE"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "MAQE");
  // pF_O
  banffnames[index] = "p_{F} O"; 
  fParLow[index] = 0.89; fParUp[index] = 1.22;
  setParName(index++, "pF_O");
  // MEC_O
  banffnames[index] = "MEC O"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "MEC_O");
  // EB_O
  banffnames[index] = "E_{B} O"; 
  fParLow[index] = 0.44; fParUp[index] = 1.56;
  setParName(index++, "EB_O");
  // CA5
  banffnames[index] = "CA5"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "CA5");
  // MA1pi
  banffnames[index] = "MANFFRES"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "MANFFRES");
  // I = 1/2 bkg
  banffnames[index] = "BgRES"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "BgRES");
  // CC other shape
  banffnames[index] = "DISMPISHP"; 
  fParLow[index] = -9999; fParUp[index] = 9999;
  setParName(index++, "DISMPISHP");
  // CC Coherent O16
  for (int i = 0; i < cccohbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("CCCOH_O %d", i);
    setParName(i+index, Form("CCCOH_O_%d", i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += cccohbins->GetNbins();
  // NC Coherent
  for (int i = 0; i < nccohbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("NCCOH %d", i);
    setParName(i+index, Form("NCCOH_%d", i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += nccohbins->GetNbins();
  // NC 1gamma
  for (int i = 0; i < nc1gammabins->GetNbins(); ++i) {
    banffnames[i+index] = Form("NC1#gamma %d", i);
    setParName(i+index, Form("NC1GAMMA_%d", i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += nc1gammabins->GetNbins();
  // NC Other
  for (int i = 0; i < ncotherbins->GetNbins(); ++i) {
    banffnames[i+index] = Form("NCOther %d", i);
    setParName(i+index, Form("NCOTHER_%d", i));
    fParLow[i+index] = 0;
    fParUp[i+index] = 9999;
  }
  index += ncotherbins->GetNbins();
  // MEC anti-nu
  banffnames[index] = "MEC #bar{#nu}"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "MEC_NUBAR");
  // CC nue/numu
  banffnames[index] = "CC#nu_{e}"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "CCNUE");
  // CC nuebar/numubar
  banffnames[index] = "CC#bar{#nu}_{e}"; 
  fParLow[index] = 0; fParUp[index] = 9999;
  setParName(index++, "CCNUEBAR");
  std::cout<<"there are "<<(size-linex[1])<<" xsec parameters."<<std::endl;

  size_xsec = size-linex[1];
  fParType = new int[size];
  fParPrior = new double[size];
  for (int i = 0; i < size; ++i) {
    if (fParName[i].find("DISMPISHP")!=std::string::npos) nominal.at(i) = 0;
    else nominal.at(i) = 1;
    fParType[i] = 0;
    fParPrior[i] = (*banff_param_prior)(i);
    fParInit[i] = fParPrior[i];
  }
  /*  
  std::cout<<"----------- nominal (unscaled) ------------\n";
  for (int i = 0; i < size; ++i) {
    std::cout << i << " " << fParName[i] << " = " << (*banff_param_nom)(i) << ", " << sqrt((*cov)(i,i)) << "\n";
  }
  std::cout<<"--------------------------------"<<std::endl;
  */
  InitPars(0.01, seed);

}

covBANFF::~covBANFF()
{
  delete fParType;
}

void covBANFF::InitPars(float stepScale, unsigned int seed)
{
  for (int i = 0; i < size; ++i) {
    fParInit[i] = (*banff_param_nom)(i);
    fParCurr[i] = fParInit[i];
    fParProp[i] = fParCurr[i];
    //fParSigma[i] = stepScale;
    fParSigma[i] = TMath::Sqrt((*cov)(i,i));
    fStepScale[i] = stepScale;
  }
  genPropKernels();
  throwNominal(true, seed);
  //randomize();
  //correlateSteps();
}

double covBANFF::getLikelihood()
{
  double LogL = 0;
  for (int i = 0; i < size; ++i) {
    if (fParProp[i] > fParUp[i] || fParProp[i] < fParLow[i]) return 99999.9;
    for (int j = 0; j < size; ++j) {
      if (fParEvalLikelihood[i] && fParEvalLikelihood[j] && fParType[i]==0 && fParType[j]==0) {
	LogL += (fParPrior[i]-fParProp[i]) * (fParPrior[i]-fParProp[i]) * (*invCov)(i,j) * 0.5;
      }
    }
  }
  return LogL;
}

double covBANFF::getPrior(int i)
{
  return fParPrior[i];
}

void covBANFF::proposeStep()
{
  randomize();
  correlateSteps();
}

void covBANFF::randomize()
{
  for (int i = 0; i < size; ++i) {
    if (fParSigma[i] > 0 && fParType[i]==0) randPar(i) = rnd->Gaus(0, 1); // fPropKernel[i]->GetRandom();
    else randPar(i) = 0;
  }
}

void covBANFF::correlateSteps()
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

void covBANFF::throwNominal(bool nomValues, unsigned int seed)
{
  TDecompChol chdcmp(*cov);
  if (!chdcmp.Decompose()) {
    std::cerr << "Cholesky decomposition failed for "<< mName << std::endl;
    exit(-1);
  }
  chel = new TMatrixD(chdcmp.GetU());
  CholeskyDecomp(size, *chel);
  ThrowParms *nom_throws = new ThrowParms(*banff_param_prior, *cov);
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
	  nominal.at(i) = (*banff_param_prior)(i);
	  continue;
	}
	if (nominal.at(i) < fParLow[i]) throw_again = true;
	if (nominal.at(i) > fParUp[i])  throw_again = true;
      }
    }
  } else {
    for (int i = 0; i < size; ++i) {
      nominal.at(i) = 1;
      if(fParName[i].find("DISMPISHP")!=std::string::npos) nominal.at(i) = 0;
    }
  }
}

double covBANFF::GetWeightFrac(int i)
{
  if (i >= 0 && i <= size) return fParProp[i] - nominal.at(i);
  else return 1;
}

double covBANFF::GetWeight(int i)
{
  if (i >= 0 && i <= size) return fParProp[i];
  else return 1;
}
