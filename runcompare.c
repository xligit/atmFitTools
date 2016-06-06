{
 gROOT->ProcessLine(".L hSplines.cxx++");
 gROOT->ProcessLine(".L histoManager.cxx++");
 gROOT->ProcessLine(".L Tool_CompareToEventByEvent.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx++");
 TChain* chmc = new TChain("h1");
 chmc->Add("./rootfiles/fake4_MC*.root");
 TTree* tr = (TTree*)chmc;
 atmFitPars* pars = new atmFitPars(3,3,7,1,"tn186");
 compareToEventByEvent* comp = new compareToEventByEvent(pars,tr,"./rootfiles/test1_factoryOutput.root","./rootfiles/test1_splineFactoryOut.root");
}
