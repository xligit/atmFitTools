{
 gROOT->ProcessLine(".L sift.C+");
 TChain ch("h1");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.07*.root");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.08*.root");
// ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.07*.root");
 sift* ss = new sift(&ch);
 ss->siftIt("postsift_atmospheric_FC_modmuoth"); 
}
