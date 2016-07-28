#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include <iostream>
#include "THStack.h"
#include "TLegend.h"
#include "hSplines.cxx"
#include "atmFitPars.h"
#include "histoTransforms.cxx"
#include "sharedPars.cxx"
#include "shared.h"
#include "splineParReader.h"
#ifndef T2K
#include "fqProcessedEvent.h"
#else
#include "t2kfqReader.h"
#include "covBase.h"
#include "covXsec.h"
#include "covBANFF.h"
#endif

//manages all histograms and splines for the fit
class histoManager{
  public:

  ///////////////////////////
  //CONSTRUCTORS//
  histoManager(int nsampl,int nbins,int ncomp,const char* name="", int nmode = 0, bool separateneutmode = false); //creates blank histogram manager
  histoManager(const char* rootfilename,int nsamp,int nbin,int ncomp,int natt, int nmode = 0, bool separateneutmode = false); //recreates a histoManager from a file
  histoManager(int nptsmc, int nptsdata); //< for unit testing, makes histoManager with gaussian histograms 
  histoManager(const char* parfile, int nmode = 0, bool separateneutmode = false); //< builds a histogram manager from histograms and values in parameter file

  ///////////////////////////
  //INTERNAL VARIABLES//
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TFile*   fin; //input file of histos to be read from memory
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  int nModes;
  TH1D* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TH1D* hMCModified[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MODIFIED MC histograms
  TH1D* hSumHisto[NSAMPMAX][NBINMAX][NATTMAX];
  TH1D* hSumHistoMod[NSAMPMAX][NBINMAX][NATTMAX];
  TH1D* hMCNeut[NSAMPMAX][NBINMAX][NCOMPMAX][NMODE][NATTMAX];
  TH1D* hMCNeutModified[NSAMPMAX][NBINMAX][NCOMPMAX][NMODE][NATTMAX];
  TH1D *hMCNeutNom[NSAMPMAX][NBINMAX][NCOMPMAX][NMODE][NATTMAX]; 
  TH1D* hMod;
  TH1D* hSum;
  TH1D* hTmp;
  int useSplineFlg;
  double normFactor;
  hSplines* theSplines[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //splines for flux/xsec params
  hSplines* moreSplines[NSAMPMAX][NBINMAX][NCOMPMAX][NMODE][NATTMAX];
  TH1D* hData[NSAMPMAX][NBINMAX][NATTMAX];  //array of all Data histograms
  TLegend* Leg;  //for histogram drawing methods
  TH2D* h2d; //for 2D debugging histograms
  double binContents[1000]; //< stores temporary bin contents for faster modificatoins
  bool separateNeutMode;
  //Array to store mean of each histogram. This affects how the histogram is scaled, since a scaling parameter
  //scales about the mean.
  double hMCMean[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];
  ///////////////////////////
  //parametrs
  atmFitPars* fitPars; 
  //void setFitPars(atmFitPars* thepars){fitPars=thepars;}
  void setFitPars(atmFitPars* thepars);
  ///////////////////////////
  //methods
  //for initialization
  void initHistos();
  void fillHistogram(int isamp, int ibin, int icomp, int iatt,double value,double weight=1.);
  void fillHistogramData(int isamp, int ibin, int iatt,double value,double weight=1.);
#ifdef T2K
  void fillHistogram(int isamp, int ibin, int icomp, int imode, int iatt, double valuem, double weight);
  void fillNominalHistogram(int isamp, int ibin, int icomp, int imode, int iatt, double value, double weight);
#endif

  ///////////////////////////
  //setters
  void setHistogram(int isamp, int ibin, int icomp, int iatt, int dataflg,TH1D* h);
#ifdef T2K
  void setHistogram(int isamp, int ibin, int icomp, int imode, int iatt, int dataflg, TH1D *h);
  void setNominalHistogram(int isamp, int ibin, int icomp, int imode, int iatt, TH1D *h);
#endif
  ///////////////////////////
  //getters
  TH1D* getHistogram(int isamp, int ibin, int icomp, int iatt);
  TH1D* getModHistogram(int isamp, int ibin, int icomp, int iatt); //gets histogram modified from atm pars
  TH1D* getHistogramData(int isamp, int ibin, int iatt){return hData[isamp][ibin][iatt];}
  hSplines* getSplines(int isamp, int ibin, int icomp, int iatt){return theSplines[isamp][ibin][icomp][iatt];}
  TH1D* getSumHistogram(int isamp, int ibin, int att, int normFlg=1);
  TH1D* getSumHistogramMod(int isamp, int ibin, int att, int normFlg=1);
  TH1D* getSplineModifiedHisto(int isamp, int ibin, int icomp, int iatt);
#ifdef T2K
  TH1D* getHistogram(int isamp, int ibin, int icomp, int imode, int iatt); // neut 
  TH1D* getNominalHistogram(int isamp, int ibin, int icomp, int imode, int iatt);
  TH1D* getModHistogram(int isamp, int ibin, int icomp, int imode, int iatt); // neut 
  hSplines* getSplines(int isamp, int ibin, int icomp, int imode, int iatt) {return moreSplines[isamp][ibin][icomp][imode][iatt];}
  TH1D* getSplineModifiedHisto(int isamp, int ibin, int icomp, int imode, int iatt);
#endif
  ///////////////////////////  
  //plotting
  void showMCBreakdown(int isample,int ibin,int iatt);
  void printBreakdownPlots(const char* directory);
  THStack* showMCBreakdownStack(int isample,int ibin,int iatt);
  void readFromFile(const char* rootename,int nsamp,int nbin,int ncomp,int natt, int nmode = 0);
  void readSplinesFromFile(const char* rootname);
  void drawSpline(int isamp, int ibin, int icomp, int iatt, int hbin, int isyst);
  void drawSpline2D( int isamp, int ibin, int icomp, int iatt, int isyst);
  ///////////////////////////
  //debugging
  void showErrorComparison(int isamp, int ibin, int iatt);
  void showSysParVariation(int isamp, int ibin, int icomp, int iatt, int isys,double varscale=1.0);
};

#endif

#ifndef HISTOMANAGER_C
#include "histoManager.cxx"
#endif
