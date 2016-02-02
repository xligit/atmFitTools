{

 ///////////////////////////////////////
 //load classes
 gROOT->ProcessLine(".L preProcess.C++");
// gROOT->ProcessLine(".L histoFactory.C++");
// gROOT->ProcessLine(".L splineFactory.C++");
 gROOT->ProcessLine(".x ~/style.c");

 preProcess* preproc = new preProcess();
 preproc->setParFileName("sharedpars.dat");
 preproc->runPreProcessing();


 ///////////////////////////////////////////////
 //name prefix for output files
// TString nameTag = "shorttest";

 ///////////////////////////////////////////
 //director for output files
 //TString directory = "./rootfiles/"; 

 ////////////////////////////////////////////////
 //setup input file chains
// TChain chmc("h1");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.00*.root");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.01*.root");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.02*.root");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.03*.root");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.04*.root");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.05*.root");
// chmc.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.06*.root");

// TChain chdata("h1");
// chdata.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.07*.root");
// chdata.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.08*.root");
// chdata.Add("/nfs/data41/t2k/amissert/atmos/fitqunRC/atmospheric/MC/feb13sk4_MassProd.reduc.09*.root");



 //////////////////////////////////////////////
 //make preProcess objects and process files

 //first for mc
// TString preProcessout = directory.Data();
// preProcessout.Append(nameTag.Data());
// preProcessout.Append("_ppmc");
// preProcess* preProcessMC = new preProcess();
 //preProcessMC->nameTag=preProcessout.Data();
 //preProcessMC->processAllFiles(&chmc);

 //now for data
 //TString preProcessout = "./rootfiles/";
 //preProcessout.Append(nameTag.Data());
 //preProcessout.Append("_ppdat");
 //preProcessMC->nameTag=preProcessout.Data();
 //preProcessMC->processAllFiles(&chdata);


}
