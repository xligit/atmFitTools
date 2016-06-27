#ifndef DEFINE_H
#define DEFINE_H

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#define T2K // should be uncommented if using T2K xsec parametrization

#define NSAMPMAX 5
#define NBINMAX 10
#define NCOMPMAX 8
#define NATTMAX 6
#ifndef T2K
#define NSYSPARMAX 20
#else
#define NSYSPARMAX 500
#define NHBINSMAX 300 // 
#define NPTSMAX 21
#define NMODE 10
#define scaling 0.1 // 0.1 for 10 years data / 100 years MC 
#endif

#endif
