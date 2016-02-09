#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include <iostream>
using namespace std;

class makeCov{
  public:
  makeCov();
  TTree* partree;
  void setParTree(TTree* tr){partree=tr;}
  TH2F* cov;
  TH2F* cor;
  void buildMatrix();
};

void makeCov::buildMatrix(){
  //setup tree stuff
  float par[1000];
  float mean[1000];
  int npar;
  partree->SetBranchAddress("par",par);
  partree->SetBranchAddress("npars",&npar);
  partree->GetEntry(0); //fills npar
  cov = new TH2F("cov","cov",npar,0.,(float)npar,npar,0.,(float)npar);
  cor = new TH2F("cor","cor",npar,0.,(float)npar,npar,0.,(float)npar);
  const int npartot = npar;
  float matrix[npartot][npartot];
  for (int i0=0;i0<npartot;i0++){
    mean[i0]=0.;
    for (int j0=0;j0<npartot;j0++){
      matrix[i0][j0] = 0.;
    }
  }
  //calc mean
  float norm=1./(float)partree->GetEntries();
 // cout<<"norm: "<<norm<<endl;
  for (int iev=0;iev<partree->GetEntries();iev++){
    partree->GetEntry(iev);
    for (int ipar=0;ipar<npartot;ipar++){
      mean[ipar]+= (par[ipar]*norm);
      //cout<<"par: "<<ipar<<" "<<par[ipar]<<endl;
    }
  }
  
  for (int kk=0;kk<npartot;kk++){
//      cout<<"mean: "<<kk<<" "<<mean[kk]<<endl;
  }

  //calc matrix
  norm = 1./((float)partree->GetEntries()-1.);
 // cout<<"norm: "<<norm<<endl;
  for (int jev=0;jev<partree->GetEntries();jev++){
    partree->GetEntry(jev);
    for (int i0=0;i0<npartot;i0++){
      for (int j0=0;j0<npartot;j0++){
        matrix[i0][j0] += ((norm)*( (par[i0]-mean[i0]) * (par[j0]-mean[j0]) ));
      }
    }
  }
  for (int kk=0;kk<npartot;kk++){
  //    cout<<"matrix: "<<kk<<" "<<matrix[kk][kk]<<endl;
  }

  //fill matrix histogram
  for (int j=0;j<npartot;j++){
    for (int k=0;k<npartot;k++){
      cov->SetBinContent(j+1,k+1,matrix[j][k]);
      cor->SetBinContent(j+1,k+1, matrix[j][k]/sqrt( (matrix[j][j]*matrix[k][k]) ));
    }
//    for (int l=0;l>j;l++){
//      cov->SetBinContent(j,l,matrix[l][j]);
//      cor->SetBinContent(j,l,matrix[l][j]/(matrix[j][j]*matrix[l][l]));
//    }
  }
  cor->Draw("colz");
  return;

}

makeCov::makeCov(){

  


}
