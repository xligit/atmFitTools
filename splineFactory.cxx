#include "splineFactory.h"

using namespace std;

const double splineFactory::sigvals[13] = {-3., -2.5, -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2., 2.5, 3.};
const double splineFactory::maqevals[21] = {-7., -6.5, -6., -5.5, -5., -4.5, -4., -3.5, -3., -2.5, -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2., 2.5, 3.};
const double splineFactory::binarys[2] = {0., 1.};

////////////////////////////////////////////////
//run spline factory using parameters from file
void splineFactory::runSplineFactory(){
  nSyst = fitPars->nSysPars;
  //create histogram manager from file with prefilled histograms
  makeManagerFromFile(runpars->hFactoryOutput.Data());
  setupHistos();
  setupSystPars();
  
  //get pointer to MC tree
  TChain chmc("h1");
  chmc.Add(runpars->hFactoryMCFiles.Data());
  setMCTree((TChain*)&chmc);
  
  //make them splines
  buildTheSplines();

  ////////////////////
  return;
}


void splineFactory::resetModHistos(){
  if (!separateNeutMode) {
    for (int ibin=0;ibin<nBin;ibin++){
      for (int isamp=0;isamp<nSamp;isamp++){
	for (int icomp=0;icomp<nComp;icomp++){
	  for (int iatt=0;iatt<nAtt;iatt++){
	    for (int ipt=0;ipt<NPTSMAX;ipt++){
	      hMC[isamp][ibin][icomp][iatt][ipt]->Reset();
	    }
	  }
	}
      }
    }
  } else {
    for (int ibin=0;ibin<nBin;ibin++){
      for (int isamp=0;isamp<nSamp;isamp++){
	for (int icomp=0;icomp<nComp;icomp++){
	  for (int iatt=0;iatt<nAtt;iatt++){
	    for (int ipt=0;ipt<NPTSMAX;ipt++){
	      hMCMode0[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode1[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode2[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode3[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode4[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode5[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode6[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode7[isamp][ibin][icomp][iatt][ipt]->Reset();
	      hMCMode8[isamp][ibin][icomp][iatt][ipt]->Reset();
	    }
	  }
	}
      }
    }
  }
}

void splineFactory::setMCTree(TChain* tr){
  mcTree = tr;
  mcEvt = new fQreader(mcTree);
  mcEvt->FillMap();
  return;
}


void splineFactory::debugtest(){
  //create new factory
  splineFactory* sfactory = new splineFactory(3,3,7,1,1,"debugtest");
  //set histogram manager and create template histograms
  sfactory->makeManagerFromFile("factoryOut_factorytest.root");
  return;
};

void splineFactory::fillLeaves(int isamp,int ibin,int icomp,int iatt,int isyst, int imode){
   //fills branches in spline information tree
   nsample = isamp;
   nbin = ibin;
   ncomponent=icomp;
   nattribute=iatt;
   nmode = imode;
   nsystpar = isyst;

   if (sysName[isyst].find("RPA_O") != std::string::npos || sysName[isyst].find("HAD") != std::string::npos) {
     npoints = 2;
     for (int i = 0; i < 2; ++i) {
       incrementSystPars(binarys[i], isyst);
       systParValues[i] = sysPar[isyst];
     }
   } else if (sysName[isyst].find("MAQE") != std::string::npos) {
     npoints = 21;
     for (int i = 0; i < 21; ++i) {
       incrementSystPars(maqevals[i], isyst);
       systParValues[i] = sysPar[isyst];
     }
   } else {
     npoints = 13;
     for (int i = 0; i < 13; ++i) {
       incrementSystPars(sigvals[i], isyst);
       systParValues[i] = sysPar[isyst];
     }
   }

   if (!separateNeutMode) {
     nhistobins = hMC[isamp][ibin][icomp][iatt][0]->GetNbinsX();
     for (int ipt=0;ipt<npoints;ipt++){
       for (int jhistobin=1;jhistobin<=nhistobins;jhistobin++){
	 if (hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin)<1e-4){
	   binWeight[ipt][jhistobin] = 1.;
	 }
	 else{
	   binWeight[ipt][jhistobin] =
	     ((double)hMC[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
	     (double)hManager->getHistogram(isamp,ibin,icomp,iatt)->GetBinContent(jhistobin);
	 }
       }
     }
   } else {
     nhistobins = hMCMode0[isamp][ibin][icomp][iatt][0]->GetNbinsX();
     for (int ipt=0;ipt<npoints;ipt++){
       for (int jhistobin=1;jhistobin<=nhistobins;jhistobin++){
	 if (hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin)<1e-4){
	   binWeight[ipt][jhistobin] = 1.;
	 } else{
	   switch (imode) {
	   case 0:
	     binWeight[ipt][jhistobin] =
	       ((double)hMCMode0[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
	       (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
	     break;
	   case 1:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode1[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
	   case 2:  
	     binWeight[ipt][jhistobin] =
               ((double)hMCMode2[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
           case 3:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode3[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
           case 4:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode4[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
           case 5:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode5[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
           case 6:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode6[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
           case 7:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode7[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
           case 8:
             binWeight[ipt][jhistobin] =
               ((double)hMCMode8[isamp][ibin][icomp][iatt][ipt]->GetBinContent(jhistobin))/
               (double)hManager->getHistogram(isamp,ibin,icomp,imode,iatt)->GetBinContent(jhistobin);
             break;
	   default:
	     binWeight[ipt][jhistobin] = 1.;
	     break;
	   }
	 }
       }
     }
   }

   /*
   double *sigvals = fitPars->getParameters();
   for (int jpt=0;jpt<NPTSMAX;jpt++){
     incrementSystPars(sigvals[jpt]);
  //   cout<<"par: "<<isyst<<sysPar[isyst]<<endl;
     systParValues[jpt]=sysPar[isyst];
   }
   return;
   */
}

void splineFactory::buildTheSplines(){

  // file setup
  TString fname = nameTag.Data();
  fname.Append("_splineOut.root");
  fout = new TFile(foutName.Data(),"recreate");

  //setup the output tree
  splineTree = new TTree("splinePars","spinePars");
  splineTree->Branch("nbin",&nbin,"nbin/I");
  splineTree->Branch("nhistobins",&nhistobins,"nhistobins/I");
  splineTree->Branch("ncomponent",&ncomponent,"ncomponent/I");
  splineTree->Branch("nattribute",&nattribute,"nattribute/I");
  splineTree->Branch("nsample",&nsample,"nsample/I");
  splineTree->Branch("nmode",&nmode,"nmode/I");
  splineTree->Branch("nsystpar",&nsystpar,"nsystpar/I");
  splineTree->Branch("nsyspartot",&nSyst,"nsyspartot/I");
  splineTree->Branch("npoints",&npoints,"npoints/I");
  splineTree->Branch("systParValues",systParValues,Form("systParValues[%d]/D",NPTSMAX));
  splineTree->Branch("binWeight",binWeight,Form("binWeight[%d][%d]/D",NPTSMAX,NHBINSMAX));  

  //setup systematic deviations
  cout<<"creating splines"<<endl; 

  for (int isyst=0;isyst<nSyst;isyst++){
    if (sysName[isyst].find("_C")!=std::string::npos) continue;
    std::cout<<"======================== "<<sysName[isyst]<<" ========================"<<std::endl;
    resetModHistos();
    for (int iev=0;iev<mcTree->GetEntries();iev++){
      //if (iev>100) break; // =========================== for defug!!!!!!
      mcTree->GetEvent(iev);  //read event
      fillAttributes();  //fill all attributes
      if (sysName[isyst].find("MAQE")!=std::string::npos) {
	for (int ipt = 0; ipt < 21; ++ipt) { // CCQE
	  incrementSystPars(maqevals[ipt], isyst);	  
	  if (iev%10000==0) std::cout<<maqevals[ipt]<<", "<<sysName[isyst]<<"="<<sysPar[isyst]<<", weight="<<eventWeight<<std::endl; 
	  getEvtWeight(isyst);
	  for (int iatt=0;iatt<nAtt;iatt++) {
	    if (!separateNeutMode) hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	    else {
	      if (mcEvt->nmode==0) hMCMode0[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==1) hMCMode1[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==2) hMCMode2[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==3) hMCMode3[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==4) hMCMode4[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==5) hMCMode5[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==6) hMCMode6[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==7) hMCMode7[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==8) hMCMode8[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	    }
	  }
	}
      } else if (sysName[isyst].find("RPA_O")!=std::string::npos || sysName[isyst].find("HAD_MUL")!=std::string::npos) { // RPA or Hadron multiplicity
	for (int ipt = 0; ipt < 2; ++ipt) { // CCQE
	  incrementSystPars(binarys[ipt], isyst);	  
	  if (iev%10000==0) std::cout<<binarys[ipt]<<", "<<sysName[isyst]<<"="<<sysPar[isyst]<<", weight="<<eventWeight<<std::endl; 
	  getEvtWeight(isyst);
	  for (int iatt=0;iatt<nAtt;iatt++) {
	    if (!separateNeutMode) hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	    else {
              if (mcEvt->nmode==0) hMCMode0[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
              else if (mcEvt->nmode==1) hMCMode1[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
              else if (mcEvt->nmode==2) hMCMode2[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==3) hMCMode3[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
              else if (mcEvt->nmode==4) hMCMode4[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==5) hMCMode5[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==6) hMCMode6[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==7) hMCMode7[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==8) hMCMode8[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	    }
	  }
	}
      } else { // all others
	for (int ipt=0;ipt<13;ipt++){ // All systematic parameters but CCQE and RPA
	  incrementSystPars(sigvals[ipt], isyst);
	  getEvtWeight(isyst); 
	  if (iev%10000==0) std::cout<<sigvals[ipt]<<", "<<sysName[isyst]<<"="<<sysPar[isyst]<<", weight="<<eventWeight<<std::endl; 
	  //fill histogram fo each attribute
	  for (int iatt=0;iatt<nAtt;iatt++){
	    if (!separateNeutMode) hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	    else {
              if (mcEvt->nmode==0) hMCMode0[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
              else if (mcEvt->nmode==1) hMCMode1[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
              else if (mcEvt->nmode==2) hMCMode2[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==3) hMCMode3[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
              else if (mcEvt->nmode==4) hMCMode4[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==5) hMCMode5[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==6) hMCMode6[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==7) hMCMode7[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	      else if (mcEvt->nmode==8) hMCMode8[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][iatt][ipt]->Fill(attribute[iatt],eventWeight);
	    }
	  } //end attribute loop
	}  //end point loop
      }
    }//end event loop

    //loop over all histograms and get spline information
    //and write to tree
    for (int kbin=0;kbin<nBin;kbin++){
      for (int kcomp=0;kcomp<nComp;kcomp++){
        for (int ksamp=0;ksamp<nSamp;ksamp++){
          for (int katt=0;katt<nAtt;katt++){
	    if (!separateNeutMode) {
	      fillLeaves(ksamp,kbin,kcomp,katt,isyst);
	      splineTree->Fill();
	    } else {
	      for (int kmode=0; kmode<nMode;++kmode) {
		fillLeaves(ksamp,kbin,kcomp,katt,isyst,kmode);
		splineTree->Fill();
	      }
	    }
          }
        }
      }
    }

  }
  splineTree->Write();
  fout->Close();
  return;
}

void splineFactory::incrementSystPars(double nsig, int isyst){
  setupSystPars();
  if (sysName[isyst].find("RPA_O") != std::string::npos) sysPar[isyst] = mcEvt->rpa[nsig];
  else if (sysName[isyst].find("HAD") != std::string::npos) sysPar[isyst] = mcEvt->rpa[nsig];
  else if (sysName[isyst].find("MAQE") != std::string::npos) sysPar[isyst] = mcEvt->maqe[nsig];
  else if (sysName[isyst].find("pF_O") != std::string::npos) sysPar[isyst] = mcEvt->pf_o[nsig];
  else if (sysName[isyst].find("EB_O") != std::string::npos) sysPar[isyst] = mcEvt->eb_o[nsig];
  else if (sysName[isyst].find("CA5") != std::string::npos) sysPar[isyst] = mcEvt->ca5[nsig];
  else if (sysName[isyst].find("MANFFRES") != std::string::npos) sysPar[isyst] = mcEvt->manffres[nsig];
  else if (sysName[isyst].find("BgRES") != std::string::npos) sysPar[isyst] = mcEvt->bgres[nsig];
  else if (sysName[isyst].find("DISMPISHP") != std::string::npos) sysPar[isyst] = mcEvt->dismpishp[nsig];
  //else if (sysName[isyst].find("RPA") != std::string::npos) sysPar[isyst] = mcEvt->rpa[nsig];
  else sysPar[isyst] = /*fitPars->sysPar[isyst]+*/sysUnc[isyst]*nsig;
  //std::cout<<nsig<<", "<<sysName[isyst]<<"="<<sysPar[isyst]<<", "<<std::endl;
}

void splineFactory::incrementSystPars(double nsig){

  //reset initial parameters
  setupSystPars();

  //change systematic parameters
  for (int isyst=0;isyst<nSyst;isyst++){
    sysPar[isyst] += sysUnc[isyst]*nsig;
    // if (sysPar[isyst]<0.) sysPar[isyst]=0.;
  }
  
  return;
}

void splineFactory::setupSystPars(){

  if (!sysParType.CompareTo("tn186")){
    nSyst=0;
    //CCQE xsec norm bin 1//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 1.0;
    nSyst++;
    //CCQE xsec norm  bin 2//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.25;
    nSyst++;
    //CCQE xsec norm bin 3//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.1;
    nSyst++;
    //CCQE xsec norm bin 4//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.05;
    nSyst++;
    //SubGeV flux norm//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.25;
    nSyst++;
    //MultiGeV flux norm//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.15;
    nSyst++;
    //CCnQE xsec norm//
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.2;
    nSyst++;
    //NC xsec norm
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.2;
    nSyst++;
    //mu/e xsec ratio
    sysPar[nSyst] = 1.0;
    sysUnc[nSyst] = 0.05;
    nSyst++;
  }
  if (!sysParType.CompareTo("t2k")) {
    nSyst = fitPars->nSysPars;
    for (int i = 0; i < nSyst; ++i) {
      sysPar[i] = fitPars->sysPar[i];
      sysUnc[i] = fitPars->sysParUnc[i];
      sysName[i] = fitPars->sysParName[i];
      //std::cout<<i<<", "<<sysName[i]<<" "<<sysPar[i]<<" +- "<<sysUnc[i]<<std::endl;
      if (sysName[i].find("RPA_O") != std::string::npos || sysName[i].find("HAD") != std::string::npos) nP[i] = 2;
      else if (sysName[i].find("MAQE") != std::string::npos) nP[i] = 21;
      else nP[i] = 13;
    }
  }
}

void splineFactory::setupHistos(){
  cout<<"called setupHistos()"<<endl;
  TH1D* htemplate;
  TString hname;
  if (!separateNeutMode) {
    for (int ipt=0;ipt<NPTSMAX;ipt++){
      for (int isamp=0;isamp<nSamp;isamp++){
	for (int ibin=0;ibin<nBin;ibin++){
	  for (int icomp=0;icomp<nComp;icomp++){
	    for (int iatt=0;iatt<nAtt;iatt++){
	      htemplate = hManager->getHistogram(isamp,ibin,icomp,iatt);
	      hname = htemplate->GetName();
	      hname.Append(Form("_%d%d%d%d%d",ipt,isamp,ibin,icomp,iatt));
	      //cout<<"setting histogram template: "<<hname.Data()<<endl;
	      hMC[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
	      hMC[isamp][ibin][icomp][iatt][ipt]->Reset();
	    }
	  }
	}
      }
    }
  } else {
    for (int ipt=0;ipt<NPTSMAX;ipt++){
      for (int isamp=0;isamp<nSamp;isamp++){
	for (int ibin=0;ibin<nBin;ibin++){
	  for (int icomp=0;icomp<nComp;icomp++){
	    for (int iatt=0;iatt<nAtt;iatt++){
	      htemplate = hManager->getHistogram(isamp,ibin,icomp,0,iatt);
	      hname = htemplate->GetName();
	      hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,0,iatt));
	      //cout<<"setting histogram template: "<<hname.Data()<<endl;
	      hMCMode0[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
	      hMCMode0[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,1,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,1,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode1[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode1[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,2,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,2,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode2[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode2[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,3,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,3,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode3[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode3[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,4,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,4,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode4[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode4[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,5,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,5,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode5[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode5[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,6,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,6,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode6[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode6[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,7,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,7,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode7[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode7[isamp][ibin][icomp][iatt][ipt]->Reset();

              htemplate = hManager->getHistogram(isamp,ibin,icomp,8,iatt);
              hname = htemplate->GetName();
              hname.Append(Form("_%d%d%d%d%d%d",ipt,isamp,ibin,icomp,8,iatt));
              //cout<<"setting histogram template: "<<hname.Data()<<endl;
              hMCMode8[isamp][ibin][icomp][iatt][ipt]=(TH1D*)htemplate->Clone(hname.Data());
              hMCMode8[isamp][ibin][icomp][iatt][ipt]->Reset();
	    }
	  }
	}
      }
    }
  }
  return;
}

void splineFactory::makeManagerFromFile(const char* fname){
  if (!separateNeutMode) {
    hManager = new histoManager(fname,nSamp,nBin,nComp,nAtt);
  } else {
    hManager = new histoManager(fname,nSamp,nBin,nComp,nAtt,9,true);
  }
  //setupHistos();
  return;
}

void splineFactory::fillAttributes(){
  attribute[0] = mcEvt->fq1rnll[0][2]-mcEvt->fq1rnll[0][1];
  attribute[1] = mcEvt->fq1rnll[1][2]-mcEvt->fq1rnll[1][1];
  return;
}

void splineFactory::fillHistograms(int ipt, int isyst){
  //fills histograms after applying systematic error parameter
  
  //reset bin contents
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jsamp=0;jsamp<nSamp;jsamp++){
      for (int jatt=0;jatt<nAtt;jatt++){
        for (int jcomp=0;jcomp<nComp;jcomp++){
          hMC[jsamp][jbin][jcomp][jatt][ipt]->Reset();
        }
      }
    }
  }

  //fill new bin contents
  for (int iev=0;iev<mcTree->GetEntries();iev++){
     mcTree->GetEvent(iev);
     fillAttributes();
     getEvtWeight(isyst);
     for (int jatt=0;jatt<nAtt;jatt++){
      hMC[mcEvt->nsample][mcEvt->nbin][mcEvt->ncomponent][jatt][ipt]->Fill(attribute[jatt],eventWeight);
    }
  }

  return;
};

double splineFactory::getEvtWeight(fQreader* mcEvt,int ipar,double value){
//  double ww = 1.;
  double ww = mcEvt->evtweight;
  int absmode = TMath::Abs(mcEvt->mode);
  double Enu     = mcEvt->pmomv[0];
  int  nutype  = TMath::Abs(mcEvt->ipnu[0]);  
 
  //simple TN186 parameterization
  if (!sysParType.CompareTo("tn186")){ 
    //CCQE norm bin1 
    if (ipar==0){
      if ((absmode==1)&&(Enu<200.)) ww*=value;
    }
    //CCQE norm bin2 
    if (ipar==1){
      if ((absmode==1)&&(Enu>200.)&&(Enu<400.)) ww*=value;
    }
    //CCQE norm bin3 
    if (ipar==2){
      if ((absmode==1)&&(Enu>400.)&&(Enu<800.)) ww*=value;
    }
    //CCQE norm bin4 
    if (ipar==3){
      if ((absmode==1)&&(Enu>800.)) ww*=value;
    }
    //SubGevFlux
    if (ipar==4){
      if (Enu<1000.) ww*=value;
    }
    //MultiGeVFlux
    if (ipar==5){
      if (Enu>1000.) ww*=value;
    }
    //CCnQE
    if (ipar==6){
      if ((absmode>1)&&(absmode<30)) ww*=value;
    }
    //NC
    if (ipar==7){
      if (absmode>=30) ww*=value;
    }
    //mu2e ratio
    if (ipar==8){
      if (nutype==14) ww*=value;
    }
  }
  if (ww<0.) ww = 0.;
  eventWeight = ww;
  return ww;

}

double splineFactory::getEvtWeight(int ipar){
//  double ww = 1.;
  double ww = mcEvt->evtweight;
  int absmode = TMath::Abs(mcEvt->mode);
  int nmode = mcEvt->nmode;
  double Enu     = mcEvt->pmomv[0];
  double evis = mcEvt->fq1rmom[0][1];
  int  nutype  = TMath::Abs(mcEvt->ipnu[0]);
//  if (ipar==0){
//   if (mcEvt->ncomponent==0) ww*=sysPar[0];
//  } 
//  if (ipar==1){
//   if (mcEvt->ncomponent==1) ww*=sysPar[1];
//  }
//  if (ipar==2){
//   if (mcEvt->ncomponent==2) ww*=sysPar[2];
//  }
   
  //CCQE norm bin1 
  if (!sysParType.CompareTo("tn186")){
    if (ipar==0){
      if ((absmode==1)&&(Enu<200.)) ww*=sysPar[0];
    }
    //CCQE norm bin2 
    if (ipar==1){
      if ((absmode==1)&&(Enu>200.)&&(Enu<400.)) ww*=sysPar[1];
    }
    //CCQE norm bin3 
    if (ipar==2){
      if ((absmode==1)&&(Enu>400.)&&(Enu<800.)) ww*=sysPar[2];
    }
    //CCQE norm bin4 
    if (ipar==3){
      if ((absmode==1)&&(Enu>800.)) ww*=sysPar[3];
    }
    //SubGevFlux
    if (ipar==4){
      if (Enu<1000.) ww*=sysPar[4];
    }
    //MultiGeVFlux
    if (ipar==5){
      if (Enu>1000.) ww*=sysPar[5];
    }
    //CCnQE
    if (ipar==6){
      if ((absmode>1)&&(absmode<30)) ww*=sysPar[6];
    }
    //NC
    if (ipar==7){
      if (absmode>=30) ww*=sysPar[7];
    }
    //mu2e ratio
    if (ipar==8){
      if (nutype==14) ww*=sysPar[8];
    }
  } else if (!sysParType.CompareTo("t2k") || !sysParType.CompareTo("banff")){
    
    if (sysName[ipar].find("MAQE")!=std::string::npos) {
      if (nmode == 0) {
	ww *= mcEvt->byEv_maqe_ccqe_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<mcEvt->byEv_maqe_ccqe_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("pF_O")!=std::string::npos) {
      if (nmode == 0) {
	ww *= mcEvt->byEv_pfo_ccqe_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<mcEvt->byEv_pfo_ccqe_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("EB_O")!=std::string::npos) {
      if (nmode == 0) {
	ww *= mcEvt->byEv_ebo_ccqe_gr->Eval(sysPar[ipar], 0, "S"); 
	std::cout<<sysName[ipar]<<" "<<mcEvt->byEv_ebo_ccqe_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("CA5")!=std::string::npos) {
      if (nmode == 1) {
	ww *= mcEvt->byEv_ca5_cc1pi_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_ca5_cc1pi_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
      else if (nmode == 4) {
	ww *= mcEvt->byEv_ca5_ncpiz_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_ca5_ncpiz_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
      else if (nmode == 5) {
	ww *= mcEvt->byEv_ca5_ncpipm_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_ca5_ncpipm_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("MANFF")!=std::string::npos) {
      if (nmode == 1) {
	ww *= mcEvt->byEv_manff_cc1pi_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_manff_cc1pi_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
      else if (nmode == 4) {
	ww *= mcEvt->byEv_manff_ncpiz_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_manff_ncpiz_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
      else if (nmode == 5) {
	ww *= mcEvt->byEv_manff_ncpipm_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_manff_ncpipm_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("BgRES")!=std::string::npos) {
      if (nmode == 1) {
	ww *= mcEvt->byEv_bgscl_cc1pi_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_bgscl_cc1pi_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
      else if (nmode == 4) {
	ww *= mcEvt->byEv_bgscl_ncpiz_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_bgscl_ncpiz_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
      else if (nmode == 5) {
	ww *= mcEvt->byEv_bgscl_ncpipm_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_bgscl_ncpipm_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("DISMPISHP")!=std::string::npos) { 
      if (nmode == 3) {
	ww *= mcEvt->byEv_dismpishp_ccoth_gr->Eval(sysPar[ipar], 0, "S");
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_dismpishp_ccoth_gr->Eval(sysPar[ipar], 0, "S")<<" | ";
      }
    }
    else if (sysName[ipar].find("FLUX_SUB")!=std::string::npos) {
      if (evis<1300) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("FLUX_MUL")!=std::string::npos) { 
      if (evis>1300) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("HAD_MUL")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("RPA_O")!=std::string::npos) {
      if (nmode == 0) {
	ww *= mcEvt->byEv_rpa_ccqe_gr->Eval(sysPar[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<mcEvt->byEv_rpa_ccqe_gr->Eval(sysPar[ipar])<<" | ";
      }
    }
    else if (sysName[ipar].find("MEC_O")!=std::string::npos) {
      if (nmode == 8) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("CCNUE_0")!=std::string::npos) {
      if (nutype == 12 && nmode <4) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("CCCOH_O_0")!=std::string::npos) {
      if (nmode == 3) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("NCOTHER_0")!=std::string::npos) {
      if (nmode == 7) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("NCCOH_0")!=std::string::npos) {
      if (nmode == 6) {
	ww *= (sysPar[ipar]+fitPars->sysParNom[ipar]);
	std::cout<<sysName[ipar]<<" "<<nmode<<" "<<sysPar[ipar]<<" | ";
      }
    }
    else if (sysName[ipar].find("FSI_INEL_LO_E")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_INEL_HI_E")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_PI_PROD")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_PI_ABS")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_CEX_LO_E")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_CEX_HI_E")!=std::string::npos) {}//do nothing for now
    else {} // all other systematic parameters don't change the weights
    
  }
  std::cout<<"----------- "<<"mode = "<<nmode<<" "<<ww<<" -----------"<<std::endl;
  
  //--------------------------------------------------------------
  /*
    if (sysName[ipar].find("MAQE_O")!=std::string::npos) {
      if (nmode == 0) ww *= mcEvt->byEv_maqe_ccqe_gr->Eval(sysPar[ipar], 0, "S");
    }
    else if (sysName[ipar].find("pF_O")!=std::string::npos) {
      if (nmode == 0) ww *= mcEvt->byEv_pfo_ccqe_gr->Eval(sysPar[ipar], 0, "S");
    }
    else if (sysName[ipar].find("EB_O")!=std::string::npos) {
      if (nmode == 0) ww *= mcEvt->byEv_ebo_ccqe_gr->Eval(sysPar[ipar], 0, "S"); 
    }
    else if (sysName[ipar].find("CA5")!=std::string::npos) {
      if (nmode == 1) ww *= mcEvt->byEv_ca5_cc1pi_gr->Eval(sysPar[ipar], 0, "S");
      else if (nmode == 4) ww *= mcEvt->byEv_ca5_ncpiz_gr->Eval(sysPar[ipar], 0, "S");
      else if (nmode == 5) ww *= mcEvt->byEv_ca5_ncpipm_gr->Eval(sysPar[ipar], 0, "S");
    }
    else if (sysName[ipar].find("MANFF")!=std::string::npos) {
      if (nmode == 1) ww *= mcEvt->byEv_manff_cc1pi_gr->Eval(sysPar[ipar], 0, "S");
      else if (nmode == 4) ww *= mcEvt->byEv_manff_ncpiz_gr->Eval(sysPar[ipar], 0, "S");
      else if (nmode == 5) ww *= mcEvt->byEv_manff_ncpipm_gr->Eval(sysPar[ipar], 0, "S");
    }
    else if (sysName[ipar].find("BgRES")!=std::string::npos) {
      if (nmode == 1) ww *= mcEvt->byEv_bgscl_cc1pi_gr->Eval(sysPar[ipar], 0, "S");
      else if (nmode == 4) ww *= mcEvt->byEv_bgscl_ncpiz_gr->Eval(sysPar[ipar], 0, "S");
      else if (nmode == 5) ww *= mcEvt->byEv_bgscl_ncpipm_gr->Eval(sysPar[ipar], 0, "S");
    }
    else if (sysName[ipar].find("DISMPISHP")!=std::string::npos) { 
      if (nmode == 3) ww *= mcEvt->byEv_dismpishp_ccoth_gr->Eval(sysPar[ipar], 0, "S");
    }
    else if (sysName[ipar].find("FLUX_SUB")!=std::string::npos) {
      if (evis<1300) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("FLUX_MUL")!=std::string::npos) { 
      if (evis>1300) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("HAD_MUL")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("RPA_O")!=std::string::npos) {
      if (nmode == 0) ww *= mcEvt->byEv_rpa_ccqe_gr->Eval(sysPar[ipar]);
    }
    else if (sysName[ipar].find("MEC_O")!=std::string::npos) {
      if (nmode == 8) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("CCNUE_0")!=std::string::npos) {
      if (nutype == 12 && nmode <4) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("CCCOH_O_0")!=std::string::npos) {
      if (nmode == 3) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("NCOTHER_0")!=std::string::npos) {
      if (nmode == 7) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("NCCOH_0")!=std::string::npos) {
      if (nmode == 6) ww *= sysPar[ipar];
    }
    else if (sysName[ipar].find("FSI_INEL_LO_E")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_INEL_HI_E")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_PI_PROD")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_PI_ABS")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_CEX_LO_E")!=std::string::npos) {}//do nothing for now
    else if (sysName[ipar].find("FSI_CEX_HI_E")!=std::string::npos) {}//do nothing for now
    else {} // all other systematic parameters don't change the weights
  }
  */
  if (ww<0.) ww = 0.;
  eventWeight = ww;
  return ww;
};

////////////////////////////////
//construct from parameter file
splineFactory::splineFactory(const char*  parfile, bool separateneutmode)
  : separateNeutMode(separateneutmode)
{
  parFileName = parfile;
  cout<<"splineFactory: Parameter File: "<<parfile<<endl;
  runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
  std::cout<<"read paramters from file"<<std::endl;
  nSamp=runpars->nSamples;
  nBin=runpars->nFVBins;
  nComp=runpars->nComponents;
  nAtt=runpars->nAttributes;
  nMode=9;
  //nSyst=runpars->nSysPars;
  nameTag=runpars->globalRootName.Data();
  foutName = runpars->splineFactoryOutput.Data();
  sysParType = runpars->sysParType;
//  foutName->Append(runpars->globalRootName.Data());
 // foutName->Append("_splines.root");
  return;
}

void splineFactory::setAtmFitPars(atmFitPars *a)
{
  fitPars = a;
}

splineFactory::splineFactory(int isamp, int ibin, int icomp, int iatt, int isyst, const char* name, int nmode, bool separateneutmode)
  : nMode(nmode)
  , separateNeutMode(separateneutmode)
{
  nameTag = name;
  nSamp = isamp;
  nBin  = ibin;
  nComp = icomp;
  nAtt  = iatt;
  nSyst = isyst;
  foutName = "splineFactoryOutput.root";
}

