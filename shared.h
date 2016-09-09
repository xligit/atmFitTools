#ifndef DEFINE_H
#define DEFINE_H

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#define T2K // should be uncommented if using T2K xsec parametrization

// uncomment this line if not using xiaoyue's skimmed atmospheric MC
//#define USE_ATM_WEIGHTS 

#define NSAMPMAX 5
#define NBINMAX 10
#define NCOMPMAX 8
#define NATTMAX 6
#define FLGDEBUG 0 // set to 1 to print out some useful things
#ifndef T2K
#define NSYSPARMAX 20
#else
#define NSYSPARMAX 500
#endif
#define NHBINSMAX 300 // 
#define NPTSMAX 21
#define NMODE 10
#define SCALING 0.1 // 0.1 for 10 years data / 100 years MC 

#endif
