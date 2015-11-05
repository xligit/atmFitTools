{
 //load necessary programs
 gROOT->ProcessLine(".L sift.C+");
 gROOT->ProcessLine(".L histoFactory.C+");
 gROOT->ProcessLine(".L splineFactory.C+");
 gROOT->ProcessLine(".x ~/style.c");

 TString nameTag = "nominalRunCheck";
 
 //setup input files
 TChain chmc("h1");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.07*.root");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.06*.root");
// chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.05*.root");
// chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.04*.root");
// chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.03*.root");
 TChain chdata("h1");
 chdata.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.09*.root");
 chdata.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.08*.root");


 //add branches to input files
// sift* siftMC = new sift(&chmc);
 TString siftName=nameTag.Data();
 siftName.Append("_MC");
// siftMC->siftIt(siftName.Data());
 TString siftOutFileMC= siftName.Data();
 siftOutFileMC.Append(".root");
 
// sift* siftData = new sift(&chdata);
 siftName=nameTag.Data();
 siftName.Append("_Data");
// siftData->siftIt(siftName.Data());
 TString siftOutFileData = siftName.Data();
 siftOutFileData.Append(".root");

 //create histograms
 TFile fmc(siftOutFileMC.Data());
 TFile fdata(siftOutFileData.Data());
 TTree* trdata = (TTree*)fdata.Get("h1");
 TTree* trmc   = (TTree*)fmc.Get("h1");
 histoFactory* hfact = new histoFactory(3,3,7,nameTag.Data()); 
 hfact->addAttribute(1);
 hfact->setDataTree(trdata);
 hfact->setMCTree(trmc);
 hfact->init();
 hfact->fillHistos();
 hfact->saveToFile();
 TString hFactoryOutput = "factoryOut_";
 hFactoryOutput.Append("nameTage.Data()");
 hFactoryOutput.Append(".root");


 //create splines
 splineFactory* s = new splineFactory(3,3,7,1,1,nameTag.Data());
 s->makeManagerFromFile(hFactoryOutput.Data());
 s->setupHistos();
 s->setupSystPars();
 s->setMCTree(mctree);
 s->buildTheSplines();
 

}
