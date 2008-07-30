void FloridaStyle()
{
  TStyle *GloStyle;
  GloStyle = gStyle; // Save the global style reference

  TStyle *FloridaStyle = new TStyle("FloridaStyle","The perfect style for plots");
  gStyle = FloridaStyle;

  // Paper size
  
  FloridaStyle->SetPaperSize(gStyle->kUSLetter);

  // Canvas

  FloridaStyle->SetCanvasColor     (  0);
  FloridaStyle->SetCanvasBorderSize( 10);
  FloridaStyle->SetCanvasBorderMode(  0);
  FloridaStyle->SetCanvasDefH      (600);
  FloridaStyle->SetCanvasDefW      (550);
  FloridaStyle->SetCanvasDefX      ( 10);
  FloridaStyle->SetCanvasDefY      ( 10);

  // Pads

  FloridaStyle->SetPadColor       (    0);
  FloridaStyle->SetPadBorderSize  (   10);
  FloridaStyle->SetPadBorderMode  (    0);
  FloridaStyle->SetPadBottomMargin(0.200);
  FloridaStyle->SetPadTopMargin   (0.080);
  FloridaStyle->SetPadLeftMargin  ( 0.18);
  FloridaStyle->SetPadRightMargin ( 0.05);
  FloridaStyle->SetPadGridX       (    0);
  FloridaStyle->SetPadGridY       (    0);
  FloridaStyle->SetPadTickX       (    0);
  FloridaStyle->SetPadTickY       (    0);

  // Frames

  FloridaStyle->SetFrameFillStyle ( 0);
  FloridaStyle->SetFrameFillColor ( 0);
  FloridaStyle->SetFrameLineColor ( 1);
  FloridaStyle->SetFrameLineStyle ( 0);
  FloridaStyle->SetFrameLineWidth ( 2);
  FloridaStyle->SetFrameBorderSize(10);
  FloridaStyle->SetFrameBorderMode( 0);

  // Histograms

  FloridaStyle->SetHistFillColor(0);
  FloridaStyle->SetHistFillStyle(1);
  FloridaStyle->SetHistLineColor(1);
  FloridaStyle->SetHistLineStyle(0);
  FloridaStyle->SetHistLineWidth(1);

  // Functions

  FloridaStyle->SetFuncColor(1);
  FloridaStyle->SetFuncStyle(0);
  FloridaStyle->SetFuncWidth(1);

  // Various

  FloridaStyle->SetTickLength (-0.015,"X");
  FloridaStyle->SetTitleSize  ( 0.050,"X");
  FloridaStyle->SetTitleOffset( 1.400,"X");
  FloridaStyle->SetLabelOffset( 0.005,"X");
  FloridaStyle->SetLabelSize  ( 0.050,"X");
  FloridaStyle->SetLabelFont  ( 42   ,"X");

  FloridaStyle->SetTickLength (-0.015,"Y");
  FloridaStyle->SetTitleSize  ( 0.050,"Y");
  FloridaStyle->SetTitleOffset( 1.200,"Y");
  FloridaStyle->SetLabelOffset( 0.015,"Y");
  FloridaStyle->SetLabelSize  ( 0.050,"Y");
  FloridaStyle->SetLabelFont  ( 42   ,"Y");

  FloridaStyle->SetStatFont       (42);
  FloridaStyle->SetTitleFont      (42);
  FloridaStyle->SetTitleBorderSize( 0);
  FloridaStyle->SetTitleFillColor (10);

  // Options

  FloridaStyle->SetOptFit (0);
  FloridaStyle->SetOptStat(0);

  return();
}
