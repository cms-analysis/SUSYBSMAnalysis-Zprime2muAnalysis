#include "RooUtils.hh"

//------------------------------------------------------------------------------
// SetRooPlot
//------------------------------------------------------------------------------
void SetRooPlot(RooPlot* frame,
		Char_t*  title,
		Char_t*  units)
{
  frame->SetTitle("");

  Double_t fitRangeBinW = frame->getFitRangeBinW();
  Char_t   xtitle[100];
  Char_t   ytitle[100];
  sprintf(xtitle,"%s [%s]",title,units);
  sprintf(ytitle,"entries / %.1f %s",fitRangeBinW,units);

  frame->SetXTitle(xtitle);
  frame->SetYTitle(ytitle);

  frame->SetLabelFont  (    42, "X");
  frame->SetLabelOffset( 0.015, "X");
  frame->SetNdivisions (   505, "X");
  frame->SetTitleFont  (    42, "X");
  frame->SetTitleOffset(  1.25, "X");
  frame->SetTitleSize  (  0.05, "X");

  frame->SetLabelFont  (    42, "Y");
  frame->SetLabelOffset( 0.020, "Y");
  frame->SetNdivisions (   505, "Y");
  frame->SetTitleFont  (    42, "Y");
  frame->SetTitleOffset(   1.9, "Y");
  frame->SetTitleSize  (  0.05, "Y");

  frame->Draw();
}
