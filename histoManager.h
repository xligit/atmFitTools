#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "fqProcessedEvent.h"
#include <iostream>
#include "THStack.h"
#include "TLegend.h"
#include "hSplines.cxx"
#include "atmFitPars.h"
#include "histoTransforms.cxx"
#include "sharedPars.cxx"
#include "shared.h"
#include "splineParReader.h"

//#define NSAMPMAX 7
//#define NCOMPMAX 20
//#define NATTMAX 20
//#define NBINMAX 10


//manages all histograms and splines for the fit
class histoManager{
  public:

  ///////////////////////////
  //CONSTRUCTORS//
  histoManager(int nsampl,int nbins,int ncomp,const char* name=""); //creates blank histogram manager
  histoManager(const char* rootname,int nsamp,int nbin,int ncomp,int natt); //recreates a histoManager from a file
  histoManager(int nptsmc, int nptsdata); //< for unit testing, makes histoManager with gaussian histograms 
  histoManager(const char* parfile); //< builds a histogram manager from histograms and values in parameter file

  ///////////////////////////
  //INTERNAL VARIABLES//
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TFile*   fin; //input file of histos to be read from memory
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  TH1D* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TH1D* hMCModified[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MODIFIED MC histograms
  TH1D* hSumHisto[NSAMPMAX][NBINMAX][NATTMAX];
  TH1D* hSumHistoMod[NSAMPMAX][NBINMAX][NATTMAX];
  TH1D* hMod;
  TH1D* hSum;
  TH1D* hTmp;
  int useSplineFlg;
  double normFactor;
  hSplines* theSplines[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //splines for flux/xsec params
  TH1D* hData[NSAMPMAX][NBINMAX][NATTMAX];  //array of all Data histograms
  TLegend* Leg;  //for histogram drawing methods
  TH2D* h2d; //for 2D debugging histograms
  double binContents[1000]; //< stores temporary bin contents for faster modificatoins
  //Array to store mean of each histogram. This affects how the histogram is scaled, since a scaling parameter
  //scales about the mean.
  double hMCMean[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];
  ///////////////////////////
  //parametrs
  atmFitPars* fitPars; 
  void setFitPars(atmFitPars* thepars){fitPars=thepars;}

  ///////////////////////////
  //methods
  //for initialization
  void initHistos();
  void fillHistogram(int isamp, int ibin, int icomp, int iatt,double value,double weight=1.);
  void fillHistogramData(int isamp, int ibin, int iatt,double value,double weight=1.);

  ///////////////////////////
  //setters
  void setHistogram(int isamp, int ibin, int icomp, int iatt, int dataflg,TH1D* h);

  ///////////////////////////
  //getters
  TH1D* getHistogram(int isamp, int ibin, int icomp, int iatt);
  TH1D* getModHistogram(int isamp, int ibin, int icomp, int iatt); //gets histogram modified from atm pars
  TH1D* getHistogramData(int isamp, int ibin, int iatt){return hData[isamp][ibin][iatt];}
  hSplines* getSplines(int isamp, int ibin, int icomp, int iatt){return theSplines[isamp][ibin][icomp][iatt];}
  TH1D* getSumHistogram(int isamp, int ibin, int att, int normFlg=1);
  TH1D* getSumHistogramMod(int isamp, int ibin, int att, int normFlg=1);
  TH1D* getSplineModifiedHisto(int isamp, int ibin, int icomp, int iatt);

  ///////////////////////////  
  //plotting
  void showMCBreakdown(int isample,int ibin,int iatt);
  THStack* showMCBreakdownStack(int isample,int ibin,int iatt);
  void readFromFile(const char* rootename,int nsamp,int nbin,int ncomp,int natt);
  void readSplinesFromFile(const char* rootname, int nsyspartot);
  
  ///////////////////////////
  //debugging
  void showErrorComparison(int isamp, int ibin, int iatt);
  void showSysParVariation(int isamp, int ibin, int icomp, int iatt, int isys,double varscale=1.0);
};

#endif

#ifndef HISTOMANAGER_C
#include "histoManager.cxx"
#endif
