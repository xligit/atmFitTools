{
gROOT->ProcessLine(".L markovTools.C++");
markovTools* mc = new markovTools(2);
float tpars[2] = {1,1};
mc->setParVar(0,1.);
mc->setParVar(1,2.);
mc->proposeStep(tpars);
}
