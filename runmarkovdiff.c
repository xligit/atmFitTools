{

gROOT->ProcessLine(".L markovDiff.cxx+");

markovDiff* md = new markovDiff("mcmctree.root",5000);

}
