{
  const char * comlOption = "+";
  gROOT->ProcessLine(".include /home/xiaoyue/atmFitTools/");
  // covariance
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/shared.h%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/ThrowParms.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/covBase.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/covXsec.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/covBANFF.cxx%s", comlOption));
  // fit parameters
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/keyread.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/sharedPars.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/atmFitPars.cxx%s", comlOption));
  // pre-Process
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/FVCalculators.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/t2kfqEvent.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/visRing.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/TH2FV.cxx%s", comlOption));
  /*
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/t2kPreProcess.cxx%s", comlOption));

  // create histograms
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/histoTransforms.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/fQreader.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/hSplines.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/splineParReader.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/histoManager.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/t2kHistoFactory.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/t2kSplineFactory.cxx%s", comlOption));
  // fit-related
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/likelihood.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/markovTools.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/histoCompare.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools/Tool_CompareToEventByEvent.cxx%s", comlOption));
  */
}
