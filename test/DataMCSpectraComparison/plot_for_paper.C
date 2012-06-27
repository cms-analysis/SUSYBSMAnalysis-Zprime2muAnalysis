// root -q -x -b plot_for_paper.C ; root -q -x -b plot_for_paper.C --muons ; root -q -x -b plot_for_paper.C --cumulative ; root -q -x -b plot_for_paper.C --muons --cumulative
// root -q -x -b plot_for_paper.C --muons ; root -q -x -b plot_for_paper.C --muons --cumulative

#include <cassert>
#include <cstdlib>
#include <string>
#include "Math/QuantFuncMathCore.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveLabel.h"
#include "TROOT.h"
#include "TStyle.h"

//#define DRAW_ZPRIME
#ifdef DRAW_ZPRIME_NO
#include "zprime.C"
#endif

using namespace std;

TGraphAsymmErrors* poisson_intervalize(TH1* h, bool zero_x, int width_normalized) {
  static const double CL = 0.6827;
  static const double alpha = (1 - CL)/2;
  static const double beta  = (1 - CL)/2;

  TGraphAsymmErrors* tgae = new TGraphAsymmErrors(h);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double nobs = h->GetBinContent(i);
    if (nobs == 0)
      continue;
    double norm_fact = 1;
    if (width_normalized) {
      norm_fact = h->GetBinWidth(i)/width_normalized;
      nobs *= norm_fact;
    }
    double lower = 0.5*ROOT::Math::chisquared_quantile_c(1-alpha, 2*nobs);
    double upper = 0.5*ROOT::Math::chisquared_quantile_c(beta, 2*(nobs+1));
    if (width_normalized) {
      nobs /= norm_fact;
      lower /= norm_fact;
      upper /= norm_fact;
    }
    if (zero_x) {
      tgae->SetPointEXlow(i-1, 0);
      tgae->SetPointEXhigh(i-1, 0);
    }
    //printf("nobs: %f -> [%f, %f] %f %f\n", nobs, lower, upper, nobs-lower, upper-nobs);
    tgae->SetPointEYlow(i-1, nobs - lower);
    tgae->SetPointEYhigh(i-1, upper - nobs);
  }

  return tgae;
}

//#define EMU

#ifdef EMU
void plot_for_paper() {
#else
void plot_for_paper2() {
#endif
  TFile fff("histos_export_emu.root");
 
  TH1* dataHist = (TH1*)fff.Get("dataHist");
  TH1* promptHist = (TH1*)fff.Get("promptHist");
  TH1* zdyHist = (TH1*)fff.Get("zdyHist");
  TH1* jetsHist = (TH1*)fff.Get("jetsHist");

  TCanvas *c1 = new TCanvas("c1", "c1",5,24,600,600);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  c1->Range(-112.5,-3.952308,1137.5,5.117294);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogy();
  c1->SetTickx();
  c1->SetTicky();
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.02);
  c1->SetFrameBorderMode(0);

  promptHist->SetTitle(";m(#mu^{+}#font[42]{e}^{-}/^{}#font[42]{e}^{+}#mu^{-}) [GeV]; Events / 20 GeV");

  promptHist->SetFillColor(2);
  promptHist->SetLineColor(2);
  jetsHist->SetFillColor(5);
  jetsHist->SetLineColor(5);

  promptHist->Draw("HIST");
  jetsHist->Draw("SAME HIST");

  TGraphAsymmErrors* dataHistPI = poisson_intervalize(dataHist, true, 0);
  dataHistPI->SetMarkerSize(0.8);
  dataHistPI->SetMarkerStyle(20);
  dataHistPI->Draw("EPZ SAME"); 

  promptHist->GetXaxis()->SetTitleSize(0.047);
  promptHist->GetXaxis()->SetTitleOffset(0.9);
  promptHist->GetYaxis()->SetTitleSize(0.047);
  promptHist->GetYaxis()->SetTitleOffset(1.2);

  promptHist->GetXaxis()->SetRangeUser(120,  900);
  promptHist->GetYaxis()->SetRangeUser(0.01, 2000);

  TLegend *leg = new TLegend(0.429, 0.668, 0.944, 0.885, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->AddEntry(dataHistPI, "DATA", "EP");
  leg->AddEntry(promptHist, "t#bar{t} + other prompt leptons", "F"); 
  leg->AddEntry(jetsHist, "jets", "F");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  TPaveLabel *pl = new TPaveLabel(0.350, 0.858, 0.810, 0.958, "CMS preliminary   #sqrt{s} = 7 TeV    #int L dt = 4.7 fb^{-1}", "brNDC");
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextSize(0.35);
  pl->Draw();

  // huge crappy hack for "EP" in TLegend::AddEntry not working
  TLine ll;
  ll.DrawLineNDC(0.494, 0.822, 0.494, 0.872);

  c1->RedrawAxis();
  system("mkdir -p plots/for_zprime_paper");
  TString fn = "MuonsElectronsOppSign";
  fn = "plots/for_zprime_paper/" + fn;
  c1->SaveAs(fn + ".pdf");
  c1->SaveAs(fn + ".root");
  c1->SaveAs(fn + ".png");
  c1->SaveAs(fn + ".gif");
  c1->SaveAs(fn + ".C");

  int k = dataHist->FindBin(60);
  int i = dataHist->FindBin(120);
  int j = dataHist->FindBin(200);
  assert(promptHist->FindBin(120) == i);
  assert(jetsHist  ->FindBin(120) == i);
  assert(promptHist->FindBin(200) == j);
  assert(jetsHist  ->FindBin(200) == j);
  printf("remember, here prompt includes Z/DY\n");
  printf("60-120 GeV:\n");
  printf("data:        %.1f\n", dataHist  ->Integral(k, i-1));
  printf("prompt+jets: %.1f\n", promptHist->Integral(k, i-1));
  printf("jets:        %.1f\n", jetsHist  ->Integral(k, i-1));
  printf("> 120 GeV:\n");
  printf("data:        %.1f\n", dataHist  ->Integral(i, 1000000));
  printf("prompt+jets: %.1f\n", promptHist->Integral(i, 1000000));
  printf("jets:        %.1f\n", jetsHist  ->Integral(i, 1000000));
  printf("> 200 GeV:\n");
  printf("data:        %.1f\n", dataHist  ->Integral(j, 1000000));
  printf("prompt+jets: %.1f\n", promptHist->Integral(j, 1000000));
  printf("jets:        %.1f\n", jetsHist  ->Integral(j, 1000000));
}

#ifndef EMU
void plot_for_paper() {
#else
void plot_for_paper2() {
#endif
  bool isCHist = false; // true for cumulative hist
  bool isElectron = true;
  int binning = 1;
  bool preliminary = false;
  bool logx = true;
  bool varbin = true;

  int argc = gApplication->Argc();
  char** argv = gApplication->Argv();
  for (int i = 0; i < argc; i++) {
    string th(argv[i]);
    if (th.find("--cumulative") != string::npos)
      isCHist = true;
    if (th.find("--muons") != string::npos)
      isElectron = false;
  }

  const char* dil_string = isElectron ? "#font[42]{ee}" : "#mu^{+}#mu^{-}";
  const char* dil_string_signed = isElectron ? "#font[42]{e}^{+}#font[42]{e}^{-}" : "#mu^{+}#mu^{-}";
  TString in_fn;
  if (isCHist)
    in_fn = isElectron ? "histos_export_heep_cumulative.root" : "histos_export_cumulative.root";
  else
    in_fn = isElectron ? "histos_export_heep.root" : "histos_export.root";
  in_fn = TString::Format("export_%s%igev/", (varbin ? "varbinround" : ""), binning) + in_fn;
  TFile fff(in_fn);
 
  TH1* dataHist = (TH1*)fff.Get("dataHist");
  TH1* zdyHist = (TH1*)fff.Get("zdyHist");
  TH1* promptHist = (TH1*)fff.Get("promptHist");
  TH1* jetsHist = (TH1*)fff.Get("jetsHist");

  TH1* zprime = 0;
#ifdef DRAW_ZPRIME
  bool draw_zprime = !isCHist;
  if (isElectron) {
    gROOT->ProcessLine(".L zprime.C");
    zprime = makeZPrime();
  }
  else
    zprime = (TH1*)fff.Get("zprime");
#else
  bool draw_zprime = false;
#endif

  TCanvas *c1 = new TCanvas("c1", "c1",5,24,600,600);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  //gStyle->SetLineWidth(2);
  c1->Range(-112.5,-3.952308,1137.5,5.117294);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogy();
  c1->SetTickx();
  c1->SetTicky();
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.04);
  c1->SetTopMargin(0.02);
  c1->SetFrameBorderMode(0);

  if (!isCHist) {
    TString binninglabel = "";
    if (binning > 0) {
      if (binning > 1)
	binninglabel = TString::Format(" / %i GeV", binning);
      else
	binninglabel = " / GeV";
    }
    zdyHist->SetTitle(TString::Format(";m(%s) [GeV]; Events%s", dil_string, binninglabel.Data()));
  }
  else
    zdyHist->SetTitle(TString::Format(";m(%s) [GeV]; Events #geq m(%s)", dil_string, dil_string));

  zdyHist->SetFillColor(7);
  zdyHist->SetLineColor(1);
  promptHist->SetFillColor(2);
  promptHist->SetLineColor(1);
  jetsHist->SetFillColor(5);
  jetsHist->SetLineColor(1);
  if (draw_zprime) {
    zprime->SetLineColor(38);
    zprime->SetLineWidth(2);
  }

  zdyHist->Draw("HIST");
  promptHist->Draw("SAME HIST");
  jetsHist->Draw("SAME HIST");
  if (draw_zprime && !isCHist) zprime->Draw("SAME HIST");

  TGraphAsymmErrors* dataHistPI = poisson_intervalize(dataHist, true, ((varbin && !isCHist) ? binning : 0));
  if (logx && isCHist) {
    // too crowded
    for (int i = 55; i < dataHistPI->GetN(); i += 2)
      dataHistPI->SetPoint(i, -5,-5);
  }
  dataHistPI->SetMarkerSize(0.6);
  dataHistPI->SetMarkerStyle(20);
  dataHistPI->Draw("EPZ SAME"); 

  if (logx) zdyHist->GetXaxis()->SetMoreLogLabels();
  if (logx) zdyHist->GetXaxis()->SetNoExponent();
  zdyHist->GetXaxis()->SetTitleSize(0.047);
  zdyHist->GetXaxis()->SetLabelOffset(0.004);
  zdyHist->GetXaxis()->SetTitleOffset(0.95);
  zdyHist->GetYaxis()->SetTitleSize(0.047);
  zdyHist->GetYaxis()->SetTitleOffset(1.2);

  if (!isCHist) {
    zdyHist->GetYaxis()->SetRangeUser(2e-5, 1.5e5);
    zdyHist->GetXaxis()->SetRangeUser(60, 2500);
  }
  else {
    zdyHist->GetYaxis()->SetRangeUser(4e-3, 3e6);
    zdyHist->GetXaxis()->SetRangeUser(60, 2500);
  }

  TLegend *leg = new TLegend(0.478, 0.667, 0.919, 0.881, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->AddEntry(dataHistPI, "DATA", "EP");
  leg->AddEntry(zdyHist, TString::Format("Z/#gamma*#rightarrow%s", dil_string_signed), "F");
  //leg->AddEntry(wjetHist, "w+jet (MC)", "F");
  leg->AddEntry(promptHist, "t#bar{t} + other prompt leptons", "F"); 
  leg->AddEntry(jetsHist, isElectron ? "jets (data)" : "jets", "F");
  if (draw_zprime && !isCHist) leg->AddEntry(zprime, TString::Format("Z'_{SSM} (1 TeV) #rightarrow %s", dil_string_signed), "F");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  TPaveLabel *pl = 0;
  if (!preliminary)
    pl = new TPaveLabel(0.44, 0.85, 0.90, 0.95, TString::Format("CMS   #sqrt{s} = 7 TeV    #int L dt = 5.%i fb^{-1}", isElectron ? 0 : 3), "brNDC");
  else
    pl = new TPaveLabel(0.39, 0.85, 0.85, 0.95, TString::Format("#splitline{CMS}{preliminary}   #sqrt{s} = 7 TeV    #int L dt = 4.%i fb^{-1}", isElectron ? 7 : 9), "brNDC");
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextSize(0.35);
  pl->Draw();

  TPaveLabel *pl2 = new TPaveLabel(0.520, 0.552, 0.979, 0.652, "#splitline{bin size}{#approx 2x detector res}", "brNDC");
  pl2->SetBorderSize(0);
  pl2->SetFillColor(0);
  pl2->SetFillStyle(0);
  pl2->SetTextSize(0.35);
  if (varbin && !isCHist)
    pl2->Draw();

  // huge crappy hack for "EP" in TLegend::AddEntry not working
  TLine ll;
  ll.DrawLineNDC(0.533, draw_zprime ? 0.811 : 0.836, 0.533, draw_zprime ? 0.859 : 0.872);

  c1->RedrawAxis();
  if (logx) c1->SetLogx(1);
  c1->SetLogy(1);
  system("mkdir -p plots/for_zprime_paper");
  TString fn;
  if (isCHist)
    fn = isElectron ? "cMassHist4600pb" : "MuonsPlusMuonsMinus_cumulative_log";
  else
    fn = isElectron ? "massHist4600pb" : "MuonsPlusMuonsMinus_log";
  fn = "plots/for_zprime_paper/" + fn;
  c1->SaveAs(fn + ".pdf");
  c1->SaveAs(fn + ".root");
  c1->SaveAs(fn + ".png");
  c1->SaveAs(fn + ".gif");
  c1->SaveAs(fn + ".C");

  /*
  for (size_t i = 0; i < dataHist->GetNbinsX()+2; ++i)
    printf("data bin %i content %.4f error %.4f\n", i, dataHist->GetBinContent(i), dataHist->GetBinError(i));
  for (size_t i = 0; i < zdyHist->GetNbinsX()+2; ++i)
    printf("zdy bin %i content %.4f error %.4f\n", i, zdyHist->GetBinContent(i), zdyHist->GetBinError(i));
  for (size_t i = 0; i < promptHist->GetNbinsX()+2; ++i)
    printf("prompt bin %i content %.4f error %.4f\n", i, promptHist->GetBinContent(i), promptHist->GetBinError(i));
  for (size_t i = 0; i < jetsHist->GetNbinsX()+2; ++i)
    printf("jets bin %i content %.4f error %.4f\n", i, jetsHist->GetBinContent(i), jetsHist->GetBinError(i));
  for (size_t i = 0; i < zprime->GetNbinsX()+2; ++i)
    printf("zprime bin %i content %.4f error %.4f\n", i, zprime->GetBinContent(i), zprime->GetBinError(i));
  printf("\n");
  */

  if (!isCHist) {
    printf("120-200 GeV:\n");
    int i = dataHist->FindBin(120);
    int j = dataHist->FindBin(200);
    assert(zdyHist  ->FindBin(120) == i);
    assert(promptHist->FindBin(120) == i);
    assert(jetsHist  ->FindBin(120) == i);
    assert(zdyHist  ->FindBin(200) == j);
    assert(promptHist->FindBin(200) == j);
    assert(jetsHist  ->FindBin(200) == j);
    printf("data:            %7.1f\n",                      dataHist  ->Integral(i, j-1));
    printf("zdy+prompt+jets: %7.1f   zdy alone:    %7.1f\n", zdyHist   ->Integral(i, j-1), zdyHist   ->Integral(i, j-1) - promptHist->Integral(i, j-1));
    printf("prompt+jets:     %7.1f   prompt alone: %7.1f\n", promptHist->Integral(i, j-1), promptHist->Integral(i, j-1) - jetsHist  ->Integral(i, j-1));
    printf("jets:            %7.1f\n",                      jetsHist  ->Integral(i, j-1));
    printf("200 GeV:\n");
    printf("data:            %7.1f\n",                      dataHist  ->Integral(j, 1000000));
    printf("zdy+prompt+jets: %7.1f   zdy alone:    %7.1f\n", zdyHist   ->Integral(j, 1000000), zdyHist   ->Integral(j, 1000000) - promptHist->Integral(j, 1000000));
    printf("prompt+jets:     %7.1f   prompt alone: %7.1f\n", promptHist->Integral(j, 1000000), promptHist->Integral(j, 1000000) - jetsHist  ->Integral(j, 1000000));
    printf("jets:            %7.1f\n",                      jetsHist  ->Integral(j, 1000000));
    if (draw_zprime) {
      int iz = zprime->FindBin(805);
      int jz = zprime->FindBin(1195);
      printf("\nzprime in 805-1195 GeV: %7.1f\n", zprime->Integral(iz, jz-1));
    }
  }
}
