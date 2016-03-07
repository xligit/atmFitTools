{
gROOT->ProcessLine(".L atmFitPars.C++");
//gROOT->ProcessLine(".L hSplines.C++");
gROOT->ProcessLine(".L histoManager.C++");
gROOT->ProcessLine(".L hSplines.C++");

//gROOT->ProcessLine(".x ~/style.c");
//gROOT->ProcessLine(".L atmFitPars.C++");

int nbin=3;
int ncomp=1;
int nsamp=1;
int natt=2;

//atmFitPars* fitpars = new atmFitpars(nbin,ncomp,nsamp,natt,0);
//TChain chdat("h1");
//TChain chmc("h1");
//chdat.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/data/patmue.run073010.004*");
//chmc.Add("/nfs/data40/t2k/amissert/skdata/cosmicsMultiGeV/14a/mc/stopmu_mc_1390_for_atmnu_test_0649*001*fQ.root");
//TChain *chdata = new TChain("h1");
//TChain *chmc   = new TChain("h1");
//chdata->Add("./rootfiles/atmfit_ppdata*.root");
//chmc->Add("./rootfiles/atmfit_ppmc*.root");
//chdata->Add("./rootfiles/nominal2_Data*.root");
//chmc->Add("./rootfiles/nominal2_MC*.root");
atmFitPars* fitpars = new atmFitPars("cosmicpars.dat");
//fitpars->initPars("tn186");
//TTree* trdata = (TTree*)chdata;
//TTree* trmc   = (TTree*)trdata;


//histoManager* hm = new histoManager(1000,10000);

//histoManager* hm = new histoManager("./rootfiles/nom3_factoryOutput.root",3,3,7,1); 
//hm->readSplinesFromFile("./rootfiles/nom2_splineFactoryOut.root",9);


//histoManager* hm = new histoManager("./rootfiles/cosmicmu_weight_histograms.root",1,3,1,2); 
//hm->setFitPars(fitpars);
//hm->readSplinesFromFile("./rootfiles/cosmic_splines.root",3);

histoManager* hm = new histoManager("cosmicpars.dat");
 

//hm->setFitPars(fitpars);

//hm->readSplinesFromFile("./rootfiles/nom2_splineFactoryOut.root",9);


//hm->setFitPars(fitpars);

}
