#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H

#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "fQreader.C"
#include <iostream>
#include "THStack.h"
#include "TLegend.h"
#include "hSplines.C"
#include "atmFitPars.C"
#include "histoTransforms.C"

#define NSAMPMAX 5
#define NCOMPMAX 20
#define NATTMAX 20
#define NBINMAX 10


//manages all histograms and splines for the fit
class histoManager{
  public:
  //CONSTRUCTORS//
  histoManager(int nsampl,int nbins,int ncomp,const char* name=""); //creates blank histogram manager
  histoManager(const char* rootname,int nsamp,int nbin,int ncomp,int natt); //recreates a histoManager from a file
  
  //INTERNAL VARIABLES//
  TString nameTag; //name associated with this instance
  TFile*   fout; //output file for filled histograms
  TFile*   fin; //input file of histos to be read from memory
  int nSamples; //number of data samples
  int nComponents; //number of MC components
  int nAttributes; //number of attributes (fiTQun outputs)
  int nBins;  //number of bins in data
  TH1F* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MC histograms
  TH1F* hMCModified[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //array of all MODIFIED MC histograms
  TH1F* hMod;
  TH1F* hSum;
  int useSplineFlg;
  hSplines* theSplines[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX]; //splines for flux/xsec params
  TH1F* hData[NSAMPMAX][NBINMAX][NATTMAX];  //array of all Data histograms
  TLegend* Leg;  //for histogram drawing methods

  //parametrs
  atmFitPars* fitPars; 
  void setFitPars(atmFitPars* thepars){fitPars=thepars;}
  //METHODS//
  void fillHistogram(int isamp, int ibin, int icomp, int iatt,float value);
  void fillHistogramData(int isamp, int ibin, int iatt,float value);

  //setters
  void setHistogram(int isamp, int ibin, int icomp, int iatt, int dataflg,TH1F* h);

  //getters
  TH1F* getHistogram(int isamp, int ibin, int icomp, int iatt);
  TH1F* getModHistogram(int isamp, int ibin, int icomp, int iatt); //gets histogram modified from atm pars
  TH1F* getHistogramData(int isamp, int ibin, int iatt){return hData[isamp][ibin][iatt];}
  hSplines* getSplines(int isamp, int ibin, int icomp, int iatt){return theSplines[isamp][ibin][icomp][iatt];}
  TH1F* getSumHistogram(int isamp, int ibin, int att);
  TH1F* getSumHistogramMod(int isamp, int ibin, int att);
  
  //plotting
  void showMCBreakdown(int isample,int ibin,int iatt);
  THStack* showMCBreakdownStack(int isample,int ibin,int iatt);
  void readFromFile(const char* rootename,int nsamp,int nbin,int ncomp,int natt);
  void readSplinesFromFile(const char* rootname, int nsyspartot);
};

#endif
