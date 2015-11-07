{
 gROOT->ProcessLine(".L sift.C+");
 TChain ch("h1");
 //ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.09*.root");
 //ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.08*.root");

 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.07*.root");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.06*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.05*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.04*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.03*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.02*.root")
 sift* ss = new sift("./rootfiles/nominal2");
 ss->processAllFiles(&ch);

// ss->processFile("/nfs/data40/t2k/amissert/skdata/atmosMC/jan14sk4_skdetsim13p90_neut532.reduc.096_fQv4r0.root");
// ss->processFile("/nfs/data40/t2k/amissert/skdata/atmosMC/jan14sk4_skdetsim13p90_neut532.reduc.084_fQv4r0.root");

// sift* ss = new sift(&ch,"nominal2");
// ss->siftIt(); 
}
