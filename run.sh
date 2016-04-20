#!/bin/bash

# run skmtuple on SKMC to calculate oscillation weight, flux weight and copy a few useful branches

# run T2KReWeight on SKMC to generate weights for xsec parameters described by splines
# This should be re-run whenever there is an update on NEUT or NIWG, or any additional spline-related xsec parameters are added in BANFF
# Either update T2KReWeight or add new dials

# run script on the aforementioned weights file to generate event-by-event splines
# This should be re-run whenever the weights files are updated
# May need to add new TGraphs in the TTree

# run macro on skimmed SKMC files and event-by-event spline files to generate final skimmed/merged files (atmFit-format files)
# This should be re-run whenever the event-by-event splilnes are updated
# May need to add new TGraphs in the TTree in preProcess; NMODE may be changed, in which case definition of nmode needs to be changed too

# run macro on atmFit-format files to generate nominal histograms and binned splines
# This should be re-run whenever BANFF or NEUT/NIWG is updated
# May need to modify covBANFF (may need to modify the sequence, change the names or set boundries mannualy ), fQreader (more TGraphs, more fQreader::FillMap()), splineFactory (splineFactory::getEvtWeight(), splineFactory::fillLeaves(), splineFacory::buildSpline())

# tune MCMC

# run fit on nominal histogram file and binned spline file

# process posterior
