{
 //load necessary programs
 gROOT->ProcessLine(".L sift.C+");
 gROOT->ProcessLine(".L histoFactory.C+");
 gROOT->ProcessLine(".L splineFactory.C+");
 gROOT->ProcessLine(".x ~/style.c");

 TString nameTag = "nominalRun";
 TString directory = "./rootfiles/"; 
 //setup input files
 TChain chmc("h1");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.07*.root");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.06*.root");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.05*.root");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.04*.root");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.03*.root");
 chmc.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.02*.root");

 TChain chdata("h1");
 chdata.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.09*.root");
 chdata.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.08*.root");


 //add branches to input files
 TString siftout = directory.Data();
 siftout.Append(nameTag.Data());
 siftout.Append("_MC");
 sift* siftMC = new sift(siftout.Data());
 siftMC->processAllFiles(&chmc);
 TString siftout = "./rootfiles/";
 siftout.Append(nameTag.Data());
 siftout.Append("_Data");
 sift* siftData = new sift(siftout.Data());
 siftData->processAllFiles(&chdata);


 

 //create histograms
 TChain *postmc=new TChain("h1");
 TChain *postdata=new TChain("h1");
 TString mcfilenames = siftMC->getFileRootName();
 mcfilenames.Append("*.root");
 postmc->Add(mcfilenames.Data());
 TString datafilenames = siftData->getFileRootName();
 datafilenames.Append("*.root");
 postdata->Add(datafilenames.Data());

 TTree* trdata = (TTree*)postdata;
 TTree* trmc   = (TTree*)postmc;
 TString factoryname = directory.Data();
 factoryname.Append(nameTag.Data());
 histoFactory* hfact = new histoFactory(3,3,7,factoryname.Data()); 
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
 TString splinename=directory.Data();
 splinename.Append(nameTag.Data());
 splineFactory* s = new splineFactory(3,3,7,1,1,splinename.Data());
 s->makeManagerFromFile(hfact->getOutputFileName().Data());
 s->setupHistos();
 s->setupSystPars();
 s->setMCTree(trmc);
 s->buildTheSplines();
 
}
