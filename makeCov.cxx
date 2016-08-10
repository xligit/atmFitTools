#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "atmFitPars.h"

#define NTOTALPARS 500

using namespace std;


/////////////////////////////////////////////////////////
//make covariance on correlation matrix from markov chain
class makeCov{

  public:

  ////////////////
  //constructor
  makeCov(const char* parfile = "");

  /////////////////
  //parameter file
  TString runparfile;

  /////////////////
  //tree of mcmc steps
  TTree* partree;

  ////////////////////////////////////////
  //set parametet tree
  void setParTree(TTree* tr){partree=tr;}

  ///////////////////////////////////////////
  //set the number of burn in steps to ignore
  int nburn;

  /////////////////////
  //matrix histograms
  TH2F* cov;
  TH2F* cor;

  //////////////////////
  //pull histograms
  TH1D* hpull;

  ///////////////////////
  //dividing lines
  TLine* lhoriz;
  TLine* lvert;

  ////////////////////////
  // number of parameters
  int nsyspar;
  int ntotpar;

  ////////////////////////
  //arrays
  double parmean[NTOTALPARS];
  double parsigma[NTOTALPARS];
  double pardefault[NTOTALPARS];
  double covarray[NTOTALPARS][NTOTALPARS];
  double corarray[NTOTALPARS][NTOTALPARS];

  /////////////////////////
  //build that matrix
  void buildMatrix();

  //////////////////////////////
  // draw nice correlation matrix
  void drawCor();

  ////////////////////////////////
  // print all 1d parameter distributions to a directory
  void printall1D(const char* dir);

};

void makeCov::printall1D(const char* dir){

  TString dirname = dir;

  TString basename = "h1D_parameter";
  TString branchname = "par";

  int npts = partree->GetEntries();

  //setup mcmc trees
  double par[500];
  int npar;
  partree->SetBranchAddress("par",par);
  partree->SetBranchAddress("npars",&npar);
  partree->GetEntry(0); //fills npar
  cout<<"Total # of parameters: "<<npar<<endl;
  cout<<"Total # of steps: "<<partree->GetEntries()<<endl;
  cout<<"Burn-in: "<<nburn<<endl;

  // canvas setup
  TCanvas* cc = new TCanvas("cc","cc",800,700);
  //partree->Draw("par[0]");
 // cc->Print("testplot.png");

  // print
  for (int ipar=0; ipar<npar; ipar++){
    TString bname = branchname.Data();
    bname.Append(Form("[%d]",ipar));
    partree->Draw(bname.Data());
    TString plotname = dirname.Data();
    plotname.Append(basename.Data());
    plotname.Append(Form("%d.png",ipar));
    cc->Print(plotname.Data());
  }

  return;      
}

void makeCov::drawCor(){
  
  cor->SetStats(0);
  cor->GetZaxis()->SetRangeUser(-1.,1.);
  cor->Draw("colz");

  double xmax = cor->GetXaxis()->GetXmax();
  double ymax = cor->GetYaxis()->GetXmax();
  double xmin = cor->GetXaxis()->GetXmin();
  double ymin = cor->GetYaxis()->GetXmin();
  double xsep = (double)ntotpar - (double)nsyspar;
  lhoriz = new TLine(0,xsep,xmax,xsep);
  lvert = new TLine(xsep,0,xsep,ymax);
  lvert->SetLineWidth(3);
  lhoriz->SetLineWidth(3);
  lhoriz->Draw("same");
  lvert->Draw("same");

  return;
}



void makeCov::buildMatrix(){

  //setup mcmc trees
  double par[500];
  int npar;
  partree->SetBranchAddress("par",par);
  partree->SetBranchAddress("npars",&npar);
  partree->GetEntry(0); //fills npar
  cout<<"Total # of parameters: "<<npar<<endl;
  cout<<"Total # of steps: "<<partree->GetEntries()<<endl;
  cout<<"Burn-in: "<<nburn<<endl;

  //create matrix templates
  cov = new TH2F("cov","cov",npar,0.,(double)npar,npar,0.,(double)npar);
  cor = new TH2F("cor","cor",npar,0.,(double)npar,npar,0.,(double)npar);

  //set initial values to zero
  const int npartot = npar; //< total number of parameters in MCMC cloud
  ntotpar = npar;
  /*
  double matrix[npartot][npartot];
  for (int i0=0;i0<npartot;i0++){
    mean[i0]=0.;
    for (int j0=0;j0<npartot;j0++){
      matrix[i0][j0] = 0.;
    }
  }
*/

  //calc means
  int npts = partree->GetEntries()-nburn;
  double norm=1./(double)npts;
  cout<<"norm: "<<norm<<endl;
  for (int iev=nburn; iev<partree->GetEntries(); iev++){
    partree->GetEntry(iev);  //< read in parameters from MCMC cloud
    for (int ipar=0; ipar<npartot; ipar++){
      parmean[ipar]+= (par[ipar]*norm);
    }
  }
  for (int kk=0;kk<npartot;kk++){
      cout<<"mean: "<<kk<<" "<<parmean[kk]<<endl;
  }

  //calc matrix
  norm = 1./((double)npts-1.);
  for (int jev=nburn;jev<partree->GetEntries();jev++){
    partree->GetEntry(jev);
    for (int i0=0;i0<npartot;i0++){
      for (int j0=0;j0<npartot;j0++){
        covarray[i0][j0] += ((norm)*( (par[i0]-parmean[i0]) * (par[j0]-parmean[j0]) ));
      }
    }
  }
  for (int kk=0;kk<npartot;kk++){
      cout<<"matrix: "<<kk<<" "<<covarray[kk][kk]<<endl;
  };

  //fill histogram of matrix values
  for (int j=0;j<npartot;j++){
    for (int k=0;k<npartot;k++){
      cov->SetBinContent(j+1,k+1,covarray[j][k]);
      cor->SetBinContent(j+1,k+1, ((covarray[j][k])/sqrt( (covarray[j][j]*covarray[k][k]) )));
      corarray[j][k] =  ((covarray[j][k])/sqrt( (covarray[j][j]*covarray[k][k]) ));
    }
    parsigma[j] = TMath::Sqrt( covarray[j][j] );
  }

  // set default parameter arrays
  atmFitPars* fitpars = new atmFitPars(runparfile.Data());
  nsyspar = fitpars->nSysPars;
  for (int ipar=0; ipar<npartot; ipar++){
    pardefault[ipar] = fitpars->getParameter(ipar);
  }

  // make pull histogram
  hpull = new TH1D("hpull","hpull",npartot,0,npartot);
  for (int ipar=0; ipar<npartot; ipar++){
    double pullvalue = (parmean[ipar] - pardefault[ipar])/parsigma[ipar];
    hpull->Fill(ipar, pullvalue);
  }
  // remove pull errors
  for (int ibin=1; ibin<hpull->GetNbinsX(); ibin++){
    hpull->SetBinError(ibin,0.);
  }
  // draw the correlation
  cor->SetContour(100);
  cor->Draw("colz");
  
  ///////////////////////
  return;

}

makeCov::makeCov(const char* parfile){

  runparfile = parfile;

 // fix initial array values to zero
 for (int i=0; i<NTOTALPARS; i++){
   parmean[i]=0.;
   parsigma[i]=0.;
   pardefault[i]=0.;
   for (int j=0; j<NTOTALPARS; j++){
      corarray[i][j] = 0.;
      covarray[i][j] = 0.;
   }
 }

}
