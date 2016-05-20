#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include <iostream>
using namespace std;


/////////////////////////////////////////////////////////
//make covariance on correlation matrix from markov chain
class makeCov{

  public:

  ////////////////
  //constructor
  makeCov();

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
  //matricies
  TH2F* cov;
  TH2F* cor;

  /////////////////////////
  //build that matrix
  void buildMatrix();
};



void makeCov::buildMatrix(){

  //setup mcmc trees
  double par[500];
  double mean[500];
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
  const int npartot = npar;
  double matrix[npartot][npartot];
  for (int i0=0;i0<npartot;i0++){
    mean[i0]=0.;
    for (int j0=0;j0<npartot;j0++){
      matrix[i0][j0] = 0.;
    }
  }

  //calc means
  int npts = partree->GetEntries()-nburn;
  double norm=1./(double)npts;
  cout<<"norm: "<<norm<<endl;
  for (int iev=nburn;iev<partree->GetEntries();iev++){
    partree->GetEntry(iev);
    for (int ipar=0;ipar<npartot;ipar++){
      mean[ipar]+= (par[ipar]*norm);
    }
  }
  for (int kk=0;kk<npartot;kk++){
      cout<<"mean: "<<kk<<" "<<mean[kk]<<endl;
  }

  //calc matrix
  norm = 1./((double)npts-1.);
  for (int jev=nburn;jev<partree->GetEntries();jev++){
    partree->GetEntry(jev);
    for (int i0=0;i0<npartot;i0++){
      for (int j0=0;j0<npartot;j0++){
        matrix[i0][j0] += ((norm)*( (par[i0]-mean[i0]) * (par[j0]-mean[j0]) ));
      }
    }
  }
  for (int kk=0;kk<npartot;kk++){
      cout<<"matrix: "<<kk<<" "<<matrix[kk][kk]<<endl;
  }

  //fill histogram of matrix values
  for (int j=0;j<npartot;j++){
    for (int k=0;k<npartot;k++){
      cov->SetBinContent(j+1,k+1,matrix[j][k]);
      cor->SetBinContent(j+1,k+1, ((matrix[j][k])/sqrt( (matrix[j][j]*matrix[k][k]) )));
    }
  }
  cor->SetContour(100);
  cor->Draw("colz");
  return;

}

makeCov::makeCov(){

  


}
