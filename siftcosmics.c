{
 gROOT->ProcessLine(".L preProcess.cxx+");
TChain chdat("h1");
TChain chmc("h1");
chdat.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run073010.004*");
chmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*001*fQ.root");
preProcess* sdat = new preProcess(&chdat);
sdat->preProcessIt("cosmicdata"); 
preProcess* smc = new preProcess(&chmc);
smc->preProcessIt("cosmicmc"); 

}
