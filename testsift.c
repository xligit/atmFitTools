{
 gROOT->ProcessLine(".L sift.C+");
 TChain ch("h1");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.09*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.08*.root");

 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.07*.root");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.06*.root");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.05*.root");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.04*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.03*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.02*.root")

 sift* ss = new sift(&ch);
 ss->siftIt("nominal4"); 
}
