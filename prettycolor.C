{
  TStyle *t2kStyle = new TStyle("T2K", "T2K colours");
  // use plain black on white colors                                                                                                                                                                        
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  //  t2kStyle->SetFillColor(0);                                                                                                                                                       
  t2kStyle->SetLegendBorderSize(1);
  // set the paper & margin sizes                                                                                                                                        
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.05);
  t2kStyle->SetPadRightMargin(0.16); //0.05                                                                 
  t2kStyle->SetPadBottomMargin(0.16);
  t2kStyle->SetPadLeftMargin(0.14);
  // use large Times-Roman fonts                                                 
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.05,"x");
  t2kStyle->SetTitleSize(0.06,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.06,"y");
  t2kStyle->SetLabelSize(0.05,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  t2kStyle->SetTitleOffset(0.96,"y");
  t2kStyle->SetPadGridX(true);
  t2kStyle->SetPadGridY(true);

  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette                                                                                                                                             

  const Int_t NRGBs = 3;
  //  const Int_t NRGBs = 5;                                                                                                                                                                                
  const Int_t NCont = 255;

  // Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };                                                                                                                                              
  // Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };                                                                                                                                              
  // Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };                                                                                                                                              
  // Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };                                                                                                                                  

  Double_t stops[NRGBs] = { 0.00, 0.5, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.00 };


  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  t2kStyle->SetNumberContours(NCont);
  gROOT->SetStyle("T2K");

  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetTickLength(0);
  gStyle->SetTickLength(0,"Y");

}
