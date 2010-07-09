#include <sstream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TFile.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TStyle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymmetryHelpers.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Functions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"

class AsymmetryParametrizer : public edm::EDFilter {
public:
  explicit AsymmetryParametrizer(const edm::ParameterSet&);
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  void endJob();

  void assemble();

  const edm::InputTag gen_particle_src;
  const bool internal_brem_on;

  const bool assemble_only;
  const std::string histos_fn;
  const std::string postscript_fn;

  TH1F* h_pt_dil;
  TH1F* h_rap_dil[2];
  TH1F* h_mass_dil[2];
  TH1F* h_phi_cs;
  TH1F* h_rap_mistag;
  TH1F* h_rap_nomistag;
  TH1F* h_cos_cs;
  TH1F* h_cos_mistag;
  TH2F* h2_rap_cos_mistag;
  TH2F* h2_rap_cos_nomistag;
  TH2F* h2_rap_cos_p;
  TH1F* h_mistag[6][3];
  TH2F* h2_mistag[4];
  TH2F* h2_pTrap;
  TH2F* h2_pTrap_mistag;
  TH2F* h2_cos_cs_vs_true;
};

AsymmetryParametrizer::AsymmetryParametrizer(const edm::ParameterSet& cfg)
  : gen_particle_src(cfg.getParameter<edm::InputTag>("gen_particle_src")),
    internal_brem_on(cfg.getParameter<bool>("internal_brem_on")),
    assemble_only(cfg.getParameter<bool>("assemble_only")),
    histos_fn(cfg.getParameter<std::string>("histos_fn")),
    postscript_fn(cfg.getParameter<std::string>("postscript_fn"))
{
  InitROOT();
  gStyle->SetOptDate(0);

  if (assemble_only) {
    assemble();
    throw cms::Exception("AssembleOnly") << "Done.\n";
  }

  asymFitManager.setConstants(cfg);
  const int nMistagBins = 20;

  edm::Service<TFileService> fs;

  h_pt_dil = fs->make<TH1F>("h_pt_dil", "pT^2 dilepton", 100, 0., 3.e6);

  h_rap_dil[0] = fs->make<TH1F>("h_rap_dil0", "abs(rap) dilepton", 50,  0. , 3.5);
  h_rap_dil[1] = fs->make<TH1F>("h_rap_dil1", "rap dilepton",      50, -3.5, 3.5);

  h_mass_dil[0] = fs->make<TH1F>("h_mass_dil0", "mass dil, full range", 100, asymFitManager.fit_win(0) - 50., asymFitManager.fit_win(1) + 50);
  h_mass_dil[1] = fs->make<TH1F>("h_mass_dil1", "mass dil, fit range",  100, asymFitManager.fit_win(0),       asymFitManager.fit_win(1));

  h_phi_cs = fs->make<TH1F>("h_phi_cs", "Collins-Soper Phi", 100, 0., TMath::Pi());

  h_rap_mistag   = fs->make<TH1F>("h_rap_mistag",   "rap dil, mistag=1", 50, 0., 3.5);
  h_rap_nomistag = fs->make<TH1F>("h_rap_nomistag", "rap dil, mistag=0", 50, 0., 3.5);

  h_cos_cs          = fs->make<TH1F>("h_cos_cs",          "cos#theta_{CS}",                             50, 0., 1.);
  h_cos_mistag      = fs->make<TH1F>("h_cos_mistag",      "cos#theta_{CS}, mistag^{1}#neq mistag^{2} ", 50, 0., 1.);

  h2_rap_cos_mistag   = fs->make<TH2F>("h2_rap_cos_mistag",   "rap vs cos dil, mistag=1", nMistagBins, 0., 1., nMistagBins, 0., 3.5);
  h2_rap_cos_nomistag = fs->make<TH2F>("h2_rap_cos_nomistag", "rap vs cos dil, mistag=0", nMistagBins, 0., 1., nMistagBins, 0., 3.5);
  h2_rap_cos_p        = fs->make<TH2F>("h2_rap_cos_p",        "rap vs cos dil",           nMistagBins, 0., 1., nMistagBins, 0., 3.5);

  h_mistag[0][0] = fs->make<TH1F>("h_mistag00", "rap all ",      50, 0., 3.5);
  h_mistag[0][1] = fs->make<TH1F>("h_mistag01", "rap mist==1",   50, 0., 3.5);

  h_mistag[1][0] = fs->make<TH1F>("h_mistag10", "pT all ",      50, 0., 20000.);
  h_mistag[1][1] = fs->make<TH1F>("h_mistag11", "pT mist==1",   50, 0., 20000.);

  h_mistag[2][0] = fs->make<TH1F>("h_mistag20", "mass all ",      50, asymFitManager.fit_win(0), asymFitManager.fit_win(1));
  h_mistag[2][1] = fs->make<TH1F>("h_mistag21", "mass mist==1",   50, asymFitManager.fit_win(0), asymFitManager.fit_win(1));

  h_mistag[3][0] = fs->make<TH1F>("h_mistag30", "cos_true all ",      50, 0., 1.);
  h_mistag[3][1] = fs->make<TH1F>("h_mistag31", "cos_true mist==1",   50, 0., 1.);

  h_mistag[4][0] = fs->make<TH1F>("h_mistag40", "dil pL all ",      50, 0., 5000.);
  h_mistag[4][1] = fs->make<TH1F>("h_mistag41", "dil pL mist==1",   50, 0., 5000.);

  h_mistag[5][0] = fs->make<TH1F>("h_mistag50", "q pL all ",      50, 0., 1000.);
  h_mistag[5][1] = fs->make<TH1F>("h_mistag51", "q pL mist==1",   50, 0., 1000.);

  h2_mistag[0] = fs->make<TH2F>("h2_mistag0", "q dil pL all", 50, 0., 5000., 50, 0., 5000.);
  h2_mistag[1] = fs->make<TH2F>("h2_mistag1", "q dil pL mist==1", 50, 0., 5000., 50, 0., 5000.);
  h2_mistag[2] = fs->make<TH2F>("h2_mistag2", "q dil pL all (reduced range)", 25, 0., 500., 25, 0., 5000.);
  h2_mistag[3] = fs->make<TH2F>("h2_mistag3", "q dil pL mist==1 (reduced range)", 25, 0., 500., 25, 0., 5000.);

  h2_pTrap        = fs->make<TH2F>("h2_pTrap",        "Rap vs pT",              nMistagBins, 0., 700., nMistagBins, 0., 3.5);
  h2_pTrap_mistag = fs->make<TH2F>("h2_pTrap_mistag", "Rap vs pT, mistag == 1", nMistagBins, 0., 700., nMistagBins, 0., 3.5);

  h2_cos_cs_vs_true = fs->make<TH2F>("h2_cos_cs_vs_true", "cos #theta_{true} vs cos #theta_{CS}", 50, -1., 1., 50, -1., 1.);
}

bool AsymmetryParametrizer::filter(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(gen_particle_src, gen_particles);

  AsymFitData data;
  if (!computeFitQuantities(*gen_particles, asymFitManager.doing_electrons(), internal_brem_on, data))
    return false; // if finding the Z'/etc failed, skip this event

  /*
  if (fake_input) {
    data.rapidity    = fake_rap[jentry];
    data.cos_true    = fake_cos_true[jentry];
    data.cos_cs      = fake_cos_cs[jentry];
    data.mistag_true = fake_mistag_true[jentry];
    data.mistag_cs = fake_mistag_cs[jentry];
  }
  */

  // Store full mass spectrum of sample while we're here.
  h_mass_dil[0]->Fill(data.mass); 

  // Fill events only in the specified mass range.
  if (data.mass < asymFitManager.fit_win(0) || data.mass > asymFitManager.fit_win(1))
    return false;

  if (data.mistag_true == 0)
    h_rap_nomistag->Fill(fabs(data.rapidity));
  else {
    h_rap_mistag->Fill(fabs(data.rapidity));
    h_mistag[0][1]->Fill(fabs(data.rapidity));
    h_mistag[1][1]->Fill(data.pT*data.pT);
    h_mistag[2][1]->Fill(data.mass);
    h_mistag[3][1]->Fill(fabs(data.cos_true));
    h_mistag[4][1]->Fill(fabs(data.pL));
    h_mistag[5][1]->Fill(fabs(data.qpL));
    h2_mistag[1]->Fill(fabs(data.qpL), fabs(data.pL));
    h2_mistag[3]->Fill(fabs(data.qpL), fabs(data.pL));
    h2_pTrap_mistag->Fill(fabs(data.pT*data.pT), fabs(data.rapidity));
  }
  h_mistag[0][0]->Fill(fabs(data.rapidity)); 
  h_mistag[1][0]->Fill(data.pT*data.pT);
  h_mistag[2][0]->Fill(data.mass);
  h_mistag[3][0]->Fill(fabs(data.cos_true));
  h_mistag[4][0]->Fill(fabs(data.pL));
  h_mistag[5][0]->Fill(fabs(data.qpL));
  h2_mistag[0]->Fill(fabs(data.qpL), fabs(data.pL));
  h2_mistag[2]->Fill(fabs(data.qpL), fabs(data.pL));
  h2_pTrap->Fill(fabs(data.pT*data.pT), fabs(data.rapidity));

  // 2D mistag probability histos for cos_cs and rapidity
  if (data.mistag_cs == 0)
    h2_rap_cos_nomistag->Fill(fabs(data.cos_cs), fabs(data.rapidity)); 
  else
    h2_rap_cos_mistag->Fill(fabs(data.cos_cs), fabs(data.rapidity)); 
  h2_rap_cos_p->Fill(fabs(data.cos_cs), fabs(data.rapidity));
	
  // These histos give small correction from fact of using CS cos 
  // instead of true
  if (data.mistag_true != data.mistag_cs)
    h_cos_mistag->Fill(fabs(data.cos_cs)); 
  h_cos_cs->Fill(fabs(data.cos_cs));

  // If you want to do convolution between CS cos and true
  if (data.mistag_cs == data.mistag_true) {
    if (data.mistag_true == 0)
      h2_cos_cs_vs_true->Fill(data.cos_true, data.cos_cs);
    else 
      h2_cos_cs_vs_true->Fill(data.cos_true, -data.cos_cs);
  }

  // Fill mass histos, one histo with all events in generated window,
  // one histo for events in reconstructed mass window, and another 
  // for events within pythia signal window.
  h_mass_dil[1]->Fill(data.mass);

  // Parametrization for pT is done with pT^2.
  h_pt_dil->Fill(data.pT*data.pT);

  // Fill histo with rapidity
  h_rap_dil[0]->Fill(fabs(data.rapidity));
  h_rap_dil[1]->Fill(data.rapidity);

  // Phi Collins-Soper is filled from 0. to Pi. (since probability
  // between Pi and 2Pi is same as from 0 to Pi). 
  if (data.phi_cs > TMath::Pi()) 
    h_phi_cs->Fill(data.phi_cs-TMath::Pi());
  else
    h_phi_cs->Fill(data.phi_cs);

  return true;
}

void AsymmetryParametrizer::endJob() {
  edm::Service<TFileService> fs;
  TH1D* h = fs->make<TH1D>("info", "", 4, 0, 4);
  h->SetBinContent(1, 1); // a way to keep track of how many jobs we were split into
  h->SetBinContent(2, double(internal_brem_on));
  h->SetBinContent(3, asymFitManager.peak_mass());
  h->SetBinContent(4, double(asymFitManager.mass_type()));
}

void AsymmetryParametrizer::assemble() {
  TFile* histos_file = new TFile(histos_fn.c_str(), "UPDATE");
  if (!histos_file->IsOpen()) {
    std::cerr << "could not open " << histos_fn << std::endl;
    return;
  }

  TDirectory* d = (TDirectory*)histos_file->Get("AsymmetryParametrizer");
  d->cd();

  h_pt_dil 	      = (TH1F*)d->Get("h_pt_dil"); 
  h_rap_dil[0] 	      = (TH1F*)d->Get("h_rap_dil0");
  h_rap_dil[1] 	      = (TH1F*)d->Get("h_rap_dil1");
  h_mass_dil[0]       = (TH1F*)d->Get("h_mass_dil0");
  h_mass_dil[1]       = (TH1F*)d->Get("h_mass_dil1");
  h_phi_cs 	      = (TH1F*)d->Get("h_phi_cs"); 
  h_rap_mistag 	      = (TH1F*)d->Get("h_rap_mistag"); 
  h_rap_nomistag      = (TH1F*)d->Get("h_rap_nomistag"); 
  h_cos_cs 	      = (TH1F*)d->Get("h_cos_cs"); 
  h_cos_mistag 	      = (TH1F*)d->Get("h_cos_mistag"); 
  h2_rap_cos_mistag   = (TH2F*)d->Get("h2_rap_cos_mistag"); 
  h2_rap_cos_nomistag = (TH2F*)d->Get("h2_rap_cos_nomistag"); 
  h2_rap_cos_p 	      = (TH2F*)d->Get("h2_rap_cos_p"); 
  h2_pTrap 	      = (TH2F*)d->Get("h2_pTrap"); 
  h2_pTrap_mistag     = (TH2F*)d->Get("h2_pTrap_mistag"); 
  h2_cos_cs_vs_true   = (TH2F*)d->Get("h2_cos_cs_vs_true"); 

  const TString titles[6] = { "rap", "pT", "mass", "cos_true", "dil pL", "q pL" };
  for (int i = 0; i < 6; ++i) {
    h_mistag[i][0] = (TH1F*)d->Get(TString::Format("h_mistag%i0", i));
    h_mistag[i][1] = (TH1F*)d->Get(TString::Format("h_mistag%i1", i));
    h_mistag[i][2] = (TH1F*)h_mistag[i][1]->Clone(TString::Format("h_mistag%i2", i));
    h_mistag[i][2]->SetTitle(titles[i] + " frac mist"); 
    h_mistag[i][2]->Sumw2(); 
  }

  h2_mistag[0] = (TH2F*)d->Get("h2_mistag0");
  h2_mistag[1] = (TH2F*)d->Get("h2_mistag1");
  h2_mistag[2] = (TH2F*)d->Get("h2_mistag2");
  h2_mistag[3] = (TH2F*)d->Get("h2_mistag3");

  TH1F* h_cos_mistag_prob = new TH1F("h_cos_mistag_prob","cos#theta_{CS} mistag frac", 50, 0., 1.);
  h_cos_mistag_prob->Sumw2();

  TH2F* h2_mistagProb = (TH2F*)h2_rap_cos_mistag->Clone("h2_mistagProb");
  h2_mistagProb->SetTitle("Mistag prob cos rap");
  TH2F* h2_pTrap_mistag_prob = (TH2F*)h2_pTrap_mistag->Clone("h2_pTrap_mistag_prob");
  h2_pTrap_mistag_prob->SetTitle("Mistag prob rap vs pT");
  TH1F* h_rap_mistag_prob = (TH1F*)h_rap_mistag->Clone("h_rap_mistag_prob");
  h_rap_mistag_prob->SetTitle("Mistag prob rap");
  TH1F* h_cos_true_mistag_prob = (TH1F*)h_mistag[3][0]->Clone("h_cos_true_mistag_prob");
  h_cos_true_mistag_prob->SetTitle("Mistag prob cos");
  TH1F* h_pL_mistag_prob = (TH1F*)h_mistag[4][0]->Clone("h_pL_mistag_prob");
  h_pL_mistag_prob->SetTitle("Mistag prob pL");
  TH2F* h2_pL_mistag_prob = (TH2F*)h2_mistag[2]->Clone("h2_pL_mistag_prob");
  h2_pL_mistag_prob->SetTitle("Mistag prob dil pL vs quark pL");

  double mistag_pars[6];
  double mass_pars[7];
  double rap_pars[5];
  double pt_pars[5];
  double phi_cs_pars[5];
  
  TH1D* info = (TH1D*)d->Get("info");
  int numJobs         = int(info->GetBinContent(1)); // if we ran originally via N jobs, the info histogram will have its stuff multiplied by N.
  bool internalBremOn = bool(info->GetBinContent(2));
  double peakMass     = info->GetBinContent(3)/numJobs;
  int massDistType    = int(info->GetBinContent(4))/numJobs;
  
  TCanvas *c1 = new TCanvas("c1", "", 0, 1, 500, 700);
  c1->SetTicks();

  TPostScript* ps = new TPostScript(postscript_fn.c_str(), 111);

  const int NUM_PAGES = 15;
  TPad *pad[NUM_PAGES];
  for (int i_page=0; i_page<NUM_PAGES; i_page++)
    pad[i_page] = new TPad("","", .05, .05, .95, .93);

  std::ostringstream page_print;
  int page = 0;

  TLatex ttl;

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  ttl.DrawLatex(.4, .95, "y of #mu^{+}#mu^{-} (mistag^{1})");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  TH1F* h_temp_rap = (TH1F*) h_rap_mistag->Clone();
  h_temp_rap->Add(h_rap_nomistag);
  h_temp_rap->Draw();
  pad[page]->cd(2); 
  h_rap_mistag->Draw();
  pad[page]->cd(3); 
  h_rap_nomistag->Draw();
  h_rap_mistag_prob->Sumw2();
  h_rap_mistag_prob->Divide(h_rap_mistag, h_temp_rap, 1., 1.);
  pad[page]->cd(4);  h_rap_mistag_prob->Draw();
  
  const double MISTAG_LIM = 3.0;

  //Fit rapidity efficiency to 0.5 + ax + bx^2
  TF1 *f_rapmis=new TF1("f_rapmis","0.5+[0]*x+[1]*x*x", 0. , MISTAG_LIM);
  f_rapmis->SetParameters(0., -0.2);
  f_rapmis->SetParNames("p1","p2");
  std::cout << "\n#### Fitting quadratic f(x=0)=0.5 to rapidity mistag prob" << std::endl;
  h_rap_mistag_prob->Fit("f_rapmis","VIR");
  double min_rap = f_rapmis->GetMinimumX(0., MISTAG_LIM);
  std::cout << "\n#### Minimum of function " << f_rapmis->GetMinimum(0., MISTAG_LIM)
       << " occurs at rap = " << min_rap << std::endl;
  mistag_pars[0] = f_rapmis->GetParameter(0);
  mistag_pars[1] = f_rapmis->GetParameter(1);
  mistag_pars[2] = min_rap;
  delete f_rapmis;
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  ttl.DrawLatex(.4, .95, 
	  "cos#theta^{*}_{CS} of #mu^{+}#mu^{-} (mistag^{1} #neq mistag^{2})");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  h_cos_mistag->Draw();
  pad[page]->cd(2); 
  h_cos_cs->Draw();
  pad[page]->cd(3); 
  h_cos_mistag_prob->Divide(h_cos_mistag, h_cos_cs, 1., 1.);
  h_cos_mistag_prob->Draw();  
  pad[page]->cd(4); 
  gPad->SetLogy(1);
  h_cos_mistag_prob->Draw();  
  TF1 *f_cosmis=new TF1("f_cosmis","[0]+[1]*exp(-[2]*x)", 0. , 1.);
  f_cosmis->SetParameters(.25, .25, 10.);
  f_cosmis->SetParNames("p0","p1","p2");
  f_cosmis->SetParLimits(0, 0., 1.);
  f_cosmis->SetParLimits(1, 0., 1.);
  f_cosmis->SetParLimits(2, 0., 50.);
  std::cout << "\n#### Fitting falling exp. to cos mistag prob" << std::endl;
  h_cos_mistag_prob->Fit("f_cosmis","VIR");
  mistag_pars[3] = f_cosmis->GetParameter(0);
  mistag_pars[4] = f_cosmis->GetParameter(1);
  mistag_pars[5] = f_cosmis->GetParameter(2);
  delete f_cosmis;
  page++;
  c1->Update();
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "y vs cos#theta^{*}_{CS} (mistag^{2})");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1);
  gPad->SetPhi(210); h2_rap_cos_mistag->Draw("lego2");
  pad[page]->cd(2);
  gPad->SetPhi(210); h2_rap_cos_nomistag->Draw("lego2");
  pad[page]->cd(3);
  gPad->SetPhi(210); h2_rap_cos_p->Draw("lego2");
  h2_mistagProb->Divide(h2_rap_cos_mistag, h2_rap_cos_p, 1., 1.);
  pad[page]->cd(4);  
  gPad->SetPhi(210); h2_mistagProb->Draw("lego2");
  page++;
  c1->Update();  
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);
  ttl.DrawLatex(.4, .95, "Mass of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  TF1* f_mass = 0;
  int nPars = 3;

  double par_norm = h_mass_dil[1]->Integral();
  double par_mean = peakMass;
  double par_fwhm = 2.*sqrt(par_mean);
  double mass_lo = h_mass_dil[1]->GetBinLowEdge(1);
  double mass_hi = h_mass_dil[1]->GetBinLowEdge(h_mass_dil[1]->GetNbinsX()+1);
  if (massDistType == AsymFitManager::MASS_EXP) {
    std::cout << "\n#### fitting falling exp to dil mass\n";
    f_mass = new TF1("f_mass", expBckg, mass_lo, mass_hi, nPars);
    f_mass->SetParNames("Norm",  "Slope", "Integral");
    f_mass->SetParameters(par_norm, -1, 1.);
    //f_mass->SetParLimits(0, 100., 1.e9);
    //f_mass->SetParLimits(1, 0., -10.);
    f_mass->FixParameter(2, 1.);
  }
  else if (massDistType == AsymFitManager::MASS_LOR) {
    std::cout << "\n#### fitting Lorentzian to dil mass\n";
    f_mass = new TF1("f_mass", Lorentzian, mass_lo, mass_hi, nPars);
    f_mass->SetParNames("Norm", "FWHM", "Mean");
    f_mass->SetParameters(par_norm, par_fwhm, par_mean);
  }
  else if (massDistType == AsymFitManager::MASS_LOREXP) {
    nPars = 6;
    std::cout << "\n#### fitting Lorentzian plus exp background to dil mass\n";
    f_mass = new TF1("f_mass", lorentzianPlusExpbckg, mass_lo, mass_hi, nPars);
    f_mass->SetParNames("NormSign", "FWHM", "Mean", "NormBckg",
			   "SlopeBckg", "IntBckg");
    f_mass->SetParameters(par_norm, par_fwhm, par_mean, 1000., -0.01, 1.);
    // set these limits
    f_mass->SetParLimits(0, 0, 1e9);
    f_mass->SetParLimits(1, 0., 1000.);
    f_mass->FixParameter(2, par_mean);
    f_mass->SetParLimits(3, 0., 1e9);
    f_mass->FixParameter(5, 1.);
  }
  else
    throw cms::Exception("AsymmetryParametrizer") << massDistType << " is not a known mass fit type!\n";

  h_mass_dil[0]->Draw();
  pad[page]->cd(2); h_mass_dil[1]->Fit("f_mass", "VIL", "", mass_lo, mass_hi);
  double mass_norm = f_mass->Integral(mass_lo, mass_hi);
  for (int i = 0; i < nPars; i++) mass_pars[i] = f_mass->GetParameter(i);
  mass_pars[nPars] = mass_norm;
  delete f_mass;
  page++;
  c1->Update();
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "Fits to y of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  std::cout << "\n#### Fitting revised thermalized cylinder model to rapidity h_rap_dil\n";
  TF1* f_rap_rtc = new TF1("f_rap_rtc", "[0]*(tanh(([1]*x)+[2]+[3])-tanh(([1]*x)-[2]+[3])+tanh(([1]*x)+[2]-[3])-tanh(([1]*x)-[2]-[3]))", 0., 4.0);
  f_rap_rtc->SetParameters(50., 1., 1., 1.);
  f_rap_rtc->SetParNames("p0", "p1", "p2", "p3");
  h_rap_dil[0]->Fit("f_rap_rtc", "VL", "", 0., 3.5);
  for (int i = 0; i < 4; i++) rap_pars[i] = f_rap_rtc->GetParameter(i);
  rap_pars[4] = 2.*f_rap_rtc->Integral(0., 3.5);
  delete f_rap_rtc;
  page++;
  c1->Update();
    
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "pT^{2} of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  gPad->SetLogy(1);
  h_pt_dil->Draw();
  h_pt_dil->SetName("2 term exponential fit"); h_pt_dil->Draw();
  std::cout << "\n#### Fitting 2 term exponential to dilepton pT^2" << std::endl;
  //TF1 *f_2exp = new TF1("f_2exp", 
  //"[0]*exp(-[1]*sqrt(x))*(x<5.e4)+[2]*exp(-[3]*sqrt(x))*(x>5.e4)", 0., 3.e6);
  TF1 *f_2exp = new TF1("f_2exp", 
       		"[0]*exp(-[1]*sqrt(x))+[2]*exp(-[3]*sqrt(x))", 0., 3.e6);
  f_2exp->SetParameters(1.3e6, 3.9e-2, 9.5e3, 1.2e-2);
  f_2exp->SetParNames("p0","p1","p2","p3");
  f_2exp->SetParLimits(0, 0., 1.e7);
  f_2exp->SetParLimits(1, 0., 1.);
  f_2exp->SetParLimits(2, 0., 1.e7);
  f_2exp->SetParLimits(3, 0., 1.);
  h_pt_dil->Fit("f_2exp", "VIL", "", 0., 3.e6);
  for (int i = 0; i < 4; i++) pt_pars[i] = f_2exp->GetParameter(i);
  pt_pars[4] = f_2exp->Integral(0., 3.e6);
  delete f_2exp;
  page++;
  c1->Update();
  
  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "#phi^{*}_{CS} of #mu^{+}#mu^{-}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(1, 2);
  pad[page]->cd(1); 
  gStyle->SetOptStat(0000);  h_phi_cs->Draw();
  h_phi_cs->SetName("x^y fit"); h_phi_cs->Draw();
  std::cout << "\n#### Fitting to Collins-Soper phi" << std::endl;
  TF1 *f_phics = new TF1("f_phics", "[0]+[1]*((x-[2])^8)", 0., 3.14);
  f_phics->SetParNames("p0","p1","p2");
  f_phics->SetParameter(0, 500.);
  f_phics->SetParameter(1, 1.);
  // set these limits
  f_phics->SetParLimits(0, 0, 1e6);
  f_phics->SetParLimits(1, 0, 1e6);
  f_phics->FixParameter(2, 1.57);
  h_phi_cs->Fit("f_phics", "VIL", "", 0., 3.14);
  if (internalBremOn) {
    for (int i = 0; i < 3; i++) phi_cs_pars[i] = f_phics->GetParameter(i);
    phi_cs_pars[3] = 8.;
    phi_cs_pars[4] = f_phics->Integral(0., 3.14);
  }
  else
    for (int i = 0; i < 5; i++) phi_cs_pars[i] = -999.;
  delete f_phics;
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  ttl.DrawLatex(.4, .95, "cos #theta_{true} vs cos #theta_{CS}");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->cd(0);
  h2_cos_cs_vs_true->Draw("lego2");
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  ttl.DrawLatex(.3, .95, "Mistag Probability for Various Quantities");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(3, 6);
  for (int i = 0; i < 6; i++) {
    h_mistag[i][2]->Divide(h_mistag[i][1], h_mistag[i][0], 1., 1., "B");
    pad[page]->cd(3*i+1);
    h_mistag[i][0]->Draw();
    pad[page]->cd(3*i+2);
    h_mistag[i][1]->Draw();
    pad[page]->cd(3*i+3);
    h_mistag[i][2]->Draw("hist E");
  }

  // Calculate histogram which can be used to obtain mistag probability
  h_cos_true_mistag_prob->Divide(h_mistag[3][1], h_mistag[3][0], 1., 1., "B");
  h_pL_mistag_prob->Divide(h_mistag[4][1], h_mistag[4][0], 1., 1., "B");
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  ttl.DrawLatex(.3, .95, "2D Mistag Probability (quark and dil pL)");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 3);
  for (int i = 0; i < 4; i++) {
    h2_mistag[i]->GetXaxis()->SetTitle("quark pL");
    h2_mistag[i]->GetYaxis()->SetTitle("dilepton pL");
    h2_mistag[i]->GetYaxis()->SetTitleOffset(1.6);
    pad[page]->cd(i+1); h2_mistag[i]->Draw();
  }
  h2_pL_mistag_prob->Sumw2();
  h2_pL_mistag_prob->Divide(h2_mistag[3], h2_mistag[2], 1., 1., "B");
  pad[page]->cd(5); 
  h2_pL_mistag_prob->GetXaxis()->SetTitle("quark pL");
  h2_pL_mistag_prob->GetXaxis()->SetTitleOffset(1.4);
  h2_pL_mistag_prob->GetYaxis()->SetTitle("dilepton pL");
  h2_pL_mistag_prob->GetYaxis()->SetTitleOffset(2.0);
  h2_pL_mistag_prob->GetZaxis()->SetTitle("mistag probability");
  h2_pL_mistag_prob->GetZaxis()->SetTitleOffset(1.5);
  gPad->SetPhi(210); h2_pL_mistag_prob->Draw("lego2");
  page++;
  c1->Update();

  ps->NewPage();
  c1->Clear();
  c1->cd(0);
  gStyle->SetOptStat(111111);
  ttl.DrawLatex(.3, .95, "2D Mistag Probability dil pT and rap");
  page_print.str(""); page_print << page + 1;
  ttl.DrawLatex(.9, .02, page_print.str().c_str());
  pad[page]->Draw();
  pad[page]->Divide(2, 2);
  pad[page]->cd(1); h2_pTrap_mistag->Draw();
  pad[page]->cd(2); h2_pTrap->Draw();
  pad[page]->cd(3);
  h2_pTrap_mistag_prob->Sumw2();
  h2_pTrap_mistag_prob->Divide(h2_pTrap_mistag, h2_pTrap, 1., 1., "B");
  h2_pTrap_mistag_prob->GetXaxis()->SetTitle("pT*pT");
  h2_pTrap_mistag_prob->GetXaxis()->SetTitleOffset(1.4);
  h2_pTrap_mistag_prob->GetYaxis()->SetTitle("Y");
  h2_pTrap_mistag_prob->GetYaxis()->SetTitleOffset(2.0);
  h2_pTrap_mistag_prob->GetZaxis()->SetTitle("mistag probability");
  h2_pTrap_mistag_prob->GetZaxis()->SetTitleOffset(1.5);
  gPad->SetPhi(210); h2_pTrap_mistag_prob->Draw("lego2");
  page++;
  c1->Update();

  // now write out all the calculated parameters and associated histos
  // to the cache file
  TArrayD arr;
  arr.Set(6, mistag_pars);
  d->WriteObject(&arr, "mistag_pars");
  arr.Set(7, mass_pars);
  d->WriteObject(&arr, "mass_pars");
  arr.Set(5, rap_pars);
  d->WriteObject(&arr, "rap_pars");
  arr.Set(5, pt_pars);
  d->WriteObject(&arr, "pt_pars");
  arr.Set(5, phi_cs_pars);
  d->WriteObject(&arr, "phi_cs_pars");
  
  for (int i = 0; i < 6; ++i)
    h_mistag[i][2]->Write();
  h2_mistagProb->Write();
  h_rap_mistag_prob->Write();
  h_cos_true_mistag_prob->Write(); 
  h_pL_mistag_prob->Write();
  h2_pL_mistag_prob->Write();
  h2_pTrap_mistag_prob->Write();
  
  histos_file->Close();
  delete histos_file;

  ps->Close(); 

  delete ps;
  delete c1;
}

DEFINE_FWK_MODULE(AsymmetryParametrizer);
