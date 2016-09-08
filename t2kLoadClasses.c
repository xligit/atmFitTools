{
  const char * comlOption = "+";
  gROOT->ProcessLine(".include /home/xiaoyue/atmFitTools_xligit/");
  // covariance
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/shared.h%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/ThrowParms.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/covBase.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/covXsec.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/covBANFF.cxx%s", comlOption));
  // fit parameters
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/keyread.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/sharedPars.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/atmFitPars.cxx%s", comlOption));
  // pre-Process
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/FVCalculators.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/t2kfqEvent.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/visRing.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/TH2FV.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/t2kPreProcess.cxx%s", comlOption));

  // create histograms
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/histoTransforms.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/masktools.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/t2kfqReader.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/hSplines.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/splineParReader.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/histoManager.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/t2kHistoFactory.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/t2kSplineFactory.cxx%s", comlOption));
  // fit-related
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/likelihood.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/markovTools.cxx%s", comlOption));
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/histoCompare.cxx%s", comlOption));
  /*
  gROOT->ProcessLine(Form(".L /home/xiaoyue/atmFitTools_xligit/Tool_CompareToEventByEvent.cxx%s", comlOption));
  */
}
