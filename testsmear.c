{
 gROOT->ProcessLine(".L histoTransforms.cxx++");
 TH1D* h = testTable(50,10,0);
 TH1D* hb1 = testBumpD(100,1,0,"t1"); 
 TH1D* hb2 = testBump(100,1,0,"t2"); 
}
