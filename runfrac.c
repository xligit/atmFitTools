{
 gROOT->ProcessLine(".L fracerr.C+");
 gROOT->ProcessLine(".x testhistomanager.c");
 
 
 
hm->fitPars->setParameter(1,33.3);
TH1D* h1 = hm->hMC[0][0][0][0];
TH1D* h2 = hm->getModHistogram(0,0,0,0);


}
