// root -q -x -b plot_for_paper.C ; root -q -x -b plot_for_paper.C --muons ; root -q -x -b plot_for_paper.C --cumulative ; root -q -x -b plot_for_paper.C --muons --cumulative

#include "Math/QuantFuncMathCore.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TROOT.h"
#include "TStyle.h"
#include "zprime.C"

TGraphAsymmErrors* poisson_intervalize(TH1* h, bool zero_x) {
  static const double CL = 0.6827;
  static const double alpha = (1 - CL)/2;
  static const double beta  = (1 - CL)/2;

  TGraphAsymmErrors* tgae = new TGraphAsymmErrors(h);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double nobs = h->GetBinContent(i);
    if (nobs == 0)
      continue;
    double lower = 0.5*ROOT::Math::chisquared_quantile_c(1-alpha, 2*nobs);
    double upper = 0.5*ROOT::Math::chisquared_quantile_c(beta, 2*(nobs+1));
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

void plot_for_paper() {
  bool isCHist = false; // true for cumulative hist
  bool isElectron = true;

  int argc = gApplication->Argc();
  char** argv = gApplication->Argv();
  for (int i = 0; i < argc; i++) {
    string th(argv[i]);
    if (th.find("--cumulative") != string::npos)
      isCHist = true;
    if (th.find("--muons") != string::npos)
      isElectron = false;
  }

  const char* dil_string = isElectron ? "ee" : "#mu^{+}#mu^{-}";
  const char* dil_string_signed = isElectron ? "e^{+}e^{-}" : "#mu^{+}#mu^{-}";
  TString in_fn;
  if (isCHist)
    in_fn = isElectron ? "cMassApprovAllFixed.root" : "histos_export_cumulative.root";
  else
    in_fn = isElectron ? "massApprovAllFixed.root" : "histos_export.root";
  TFile fff(in_fn);
 
  TH1* dataHist = (TH1*)fff.Get("dataHist");
  TH1* zeeHist = (TH1*)fff.Get("zeeHist");
  TH1* ttbarHist = (TH1*)fff.Get("ttbarHist");
  TH1* qcdHist = (TH1*)fff.Get("qcdHist");
  TH1* zprime = 0;

  if (isElectron) {
    gROOT->ProcessLine(".L zprime.C");
    zprime = makeZPrime();
  }
  else
    zprime = (TH1*)fff.Get("zprime");

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
  c1->SetRightMargin(0.07);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);

  if (!isCHist) {
    zeeHist->SetTitle(TString::Format(";m(%s) [GeV]; Events / 5 GeV", dil_string));
    zeeHist->GetXaxis()->SetNdivisions(505);
  }
  else
    zeeHist->SetTitle(TString::Format(";m(%s) [GeV]; Events #geq m(%s)", dil_string, dil_string));

  zeeHist->SetFillColor(7);
  zeeHist->SetLineColor(7);
  ttbarHist->SetFillColor(2);
  ttbarHist->SetLineColor(2);
  qcdHist->SetFillColor(4);
  qcdHist->SetLineColor(4);
  zprime->SetLineColor(38);
  zprime->SetLineWidth(2);

  zeeHist->Draw("HIST");
  ttbarHist->Draw("SAME HIST");
  qcdHist->Draw("SAME HIST");
  if (!isCHist) zprime->Draw("SAME HIST");

  TGraphAsymmErrors* dataHistPI = poisson_intervalize(dataHist, true);
  dataHistPI->SetMarkerSize(0.8);
  dataHistPI->SetMarkerStyle(20);
  dataHistPI->Draw("EPZ SAME"); 

  zeeHist->GetXaxis()->SetTitleSize(0.047);
  zeeHist->GetXaxis()->SetTitleOffset(0.9);
  zeeHist->GetYaxis()->SetTitleSize(0.047);
  zeeHist->GetYaxis()->SetTitleOffset(1.2);

  if (!isCHist) {
    zeeHist->GetYaxis()->SetRangeUser(1e-3, 2e4);
    zeeHist->GetXaxis()->SetRangeUser(50, 1000);
  }
  else {
    zeeHist->GetYaxis()->SetRangeUser(5e-1, 3e4);
    zeeHist->GetXaxis()->SetRangeUser(50, 500);
  }

  TLegend *leg = new TLegend(0.48, 0.53, 0.88, 0.88, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->AddEntry(dataHistPI, "DATA", "EP");
  leg->AddEntry(zeeHist, TString::Format("Z/#gamma*#rightarrow%s", dil_string_signed), "F");
  //leg->AddEntry(wjetHist, "w+jet (MC)", "F");
  leg->AddEntry(ttbarHist, "t#bar{t} + other prompt leptons", "F"); 
  leg->AddEntry(qcdHist, isElectron ? "jets (data)" : "jets", "F");
  if (!isCHist) leg->AddEntry(zprime, TString::Format("Z'_{SSM} (750 GeV) #rightarrow %s", dil_string_signed), "F");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  TPaveLabel *pl = new TPaveLabel(0.40, 0.89, 0.86, 0.99, TString::Format("CMS    #sqrt{s} = 7 TeV    #int L dt = %u pb^{-1}", isElectron ? 35 : 40), "brNDC");
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextSize(0.35);
  pl->Draw();

  // huge crappy hack for "EP" in TLegend::AddEntry not working
  TLine ll;
  if (isCHist)
    ll.DrawLineNDC(0.53, 0.805, 0.53, 0.865);
  else
    ll.DrawLineNDC(0.53, 0.815, 0.53, 0.87);

  c1->RedrawAxis();
  c1->SetLogy(1);
  system("mkdir -p plots/for_paper");
  TString fn;
  if (isCHist)
    fn = isElectron ? "cMassHist35ForPaper" : "MuonsPlusMuonsMinus_cumulative_log";
  else
    fn = isElectron ? "massHist35ForPaper" : "MuonsPlusMuonsMinus_log";
  fn = "plots/for_paper/" + fn;
  c1->SaveAs(fn + ".pdf");
  c1->SaveAs(fn + ".root");
  c1->SaveAs(fn + ".png");
  c1->SaveAs(fn + ".C");

  /*
  for (size_t i = 0; i < dataHist->GetNbinsX()+2; ++i)
    printf("data bin %i content %.4f error %.4f\n", i, dataHist->GetBinContent(i), dataHist->GetBinError(i));
  for (size_t i = 0; i < zeeHist->GetNbinsX()+2; ++i)
    printf("zee bin %i content %.4f error %.4f\n", i, zeeHist->GetBinContent(i), zeeHist->GetBinError(i));
  for (size_t i = 0; i < ttbarHist->GetNbinsX()+2; ++i)
    printf("ttbar bin %i content %.4f error %.4f\n", i, ttbarHist->GetBinContent(i), ttbarHist->GetBinError(i));
  for (size_t i = 0; i < qcdHist->GetNbinsX()+2; ++i)
    printf("qcd bin %i content %.4f error %.4f\n", i, qcdHist->GetBinContent(i), qcdHist->GetBinError(i));
  for (size_t i = 0; i < zprime->GetNbinsX()+2; ++i)
    printf("zprime bin %i content %.4f error %.4f\n", i, zprime->GetBinContent(i), zprime->GetBinError(i));
  printf("\n");
  */

  if (!isCHist) {
    printf("120-200 GeV:\n");
    int i = dataHist->FindBin(120);
    int j = dataHist->FindBin(200);
    int iz = zprime->FindBin(120);
    int jz = zprime->FindBin(200);
    assert(zeeHist  ->FindBin(120) == i);
    assert(ttbarHist->FindBin(120) == i);
    assert(qcdHist  ->FindBin(120) == i);
    assert(zprime   ->FindBin(120) == iz);
    assert(zeeHist  ->FindBin(200) == j);
    assert(ttbarHist->FindBin(200) == j);
    assert(qcdHist  ->FindBin(200) == j);
    assert(zprime   ->FindBin(200) == jz);
    printf("data:   %.1f\n", dataHist ->Integral(i, j-1));
    printf("zee:    %.1f\n", zeeHist  ->Integral(i, j-1));
    printf("ttbar:  %.1f\n", ttbarHist->Integral(i, j-1));
    printf("qcd:    %.1f\n", qcdHist  ->Integral(i, j-1));
    printf("zprime: %.1f\n", zprime   ->Integral(iz, jz-1));
    printf("200 GeV:\n");
    printf("data:   %.1f\n", dataHist ->Integral(j, 1000000));
    printf("zee:    %.1f\n", zeeHist  ->Integral(j, 1000000));
    printf("ttbar:  %.1f\n", ttbarHist->Integral(j, 1000000));
    printf("qcd:    %.1f\n", qcdHist  ->Integral(j, 1000000));
    printf("zprime: %.1f\n", zprime   ->Integral(jz, 1000000));
  }
}
