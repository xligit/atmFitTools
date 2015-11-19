{
 gROOT->ProcessLine(".L hSplines.C++");
 gROOT->ProcessLine(".L histoManager.C++");
 gROOT->ProcessLine(".L Tool_CompareToEventByEvent.C++");
 gROOT->ProcessLine(".L atmFitPars.C++");
 TChain* chmc = new TChain("h1");
 chmc->Add("./rootfiles/fake3_MC*.root");
 TTree* tr = (TTree*)chmc;
 atmFitPars* pars = new atmFitPars(3,3,7,1,"tn186");
 compareToEventByEvent* comp = new compareToEventByEvent(pars,tr,"./rootfiles/fake3_hFactoryOutput.root","./rootfiles/fake3_splineOut.root");
}
