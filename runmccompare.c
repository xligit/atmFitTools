{
  gROOT->ProcessLine(".L markovCompare.C+");

//  markovCompare* mc = new markovCompare("./mcmc/fake1_biasonly_minfit2.root","./mcmc/fake1_biasonly_truefit2.root");
//  mc->nburn1 = 50000;
//  mc->nburn2 = 50000;

//  markovCompare* mc = new markovCompare("./mcmc/fake1_allpars_mcmc.root","./mcmc/fake1_allpars_allpars_demcmc.root");
//  mc->nburn1 = 50000;
//  mc->nburn2 = 10000;

//  markovCompare* mc = new markovCompare("./mcmc/closure_allpar_mcmcfit.root","mcmctree.root");
//  mc->nburn1 = 90000;
//  mc->nburn2 = 5000;

//  markovCompare* mc = new markovCompare("./mcmctree.root","./mcmc/tgraph_biasonly.root");
  markovCompare* mc = new markovCompare("./mcmc/results/fake1_allpar.root","./mcmc/results/fake4_bin6fit.root");
  mc->nburn1 = 5000;
  mc->nburn2 = 50000000;

//  markovCompare* mc = new markovCompare("./mcmc/hpi0_mcmcfit.root","./mcmc/hpi0_demcmcfit.root");
//  mc->nburn1 = 50000;
//  mc->nburn2 = 50000;

//  markovCompare* mc = new markovCompare("./mcmc/simple1_minfit.root","./mcmc/simple1_truefit.root");
//  mc->nburn1 = 10000;
//  mc->nburn2 = 10000;
}
