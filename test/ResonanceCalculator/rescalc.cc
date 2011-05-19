#include "TFile.h"
#include "TH1D.h"
#include "PhysicsTools/RooStatsCms/interface/ResonanceCalculators.hh"

int main() {
  TFile *file = TFile::Open("ana_datamc_data.root");
  file->cd();
  TH1* hist = (TH1*)file->GetDirectory("OurMuonsPlusMuonsMinusHistos")->Get("DileptonMass");
  hist->Rebin(5);
  TH1D* newhist = new TH1D("newhist", "new hist", 100, 200, 700);
  for (int i = 1; i <= 100; ++i) {
    newhist->SetBinContent(i, hist->GetBinContent(i+40));
    newhist->SetBinError  (i, hist->GetBinError  (i+40));
  }

  FactoryResCalc rc("signal",      "RooGaussian::signal(obs, signalmass, signalwidth)",
		    "background",  "EXPR::background('exp(-0.006912*obs)/pow(obs, 2.404)', obs)",
		    "signalwidth", "expr::signalwidth('0.138e-1*signalmass + 0.9315e-4*signalmass*signalmass - 0.1077e-7*signalmass*signalmass*signalmass', signalmass)");
  rc.setBinnedData(newhist);
  rc.setNumBinsToDraw(newhist->GetNbinsX());
  rc.setMinMaxSignalMass(200., 700.);
  rc.setNumPseudoExperiments(0);
  rc.setRandomSeed(1);
  runResCalc(rc, "output");
}
