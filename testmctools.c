{
gROOT->ProcessLine(".L markovTools.cxx++");
markovTools* mc = new markovTools(2);
float tpars[2] = {1,1};
mc->setParVar(0,1.);
mc->setParVar(1,2.);
mc->proposeStep(tpars);
}
