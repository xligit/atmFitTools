{
TFile fout("testspline.root","recreate");

double x[4] = {1,2,3,4};
double y[4] = {-2,0,4,1};

TSpline3 *ts = new TSpline3("ts",x,y,4);

ts->SetName("testspline");
fout.Write();

}
