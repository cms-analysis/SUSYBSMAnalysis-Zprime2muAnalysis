#include <iomanip>

#include "TFile.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFitManager.h"

void AsymFitManager::setConstants(const edm::ParameterSet& pset) {
  if (_mistag_calc) {
    delete _mistag_calc;
    _mistag_calc = 0;
  }

  _doing_electrons = pset.getParameter<bool>("doingElectrons");
  _lepton_mass = _doing_electrons ? 0.000511 : 0.10566;

  _mass_type = pset.getParameter<int>("massDistType");
  if (_mass_type != MASS_EXP && _mass_type != MASS_LOR && _mass_type != MASS_LOREXP)
    throw cms::Exception("BadMassDistType") << "massDistType = " << _mass_type << "\n";

  _max_rap = pset.getParameter<double>("maxRapidity");
  _max_pt = pset.getParameter<double>("maxPt");
  _fit_win = pset.getParameter<std::vector<double> >("fitWindow");
  _peak_mass = pset.getParameter<double>("peakMass");
  _rec_sigma = pset.getParameter<std::vector<double> >("recSigma");

  _correct_mistags  = pset.getParameter<bool>("correctMistags");
  _calculate_mistag = pset.getParameter<bool>("calculateMistag");
  _use_mistag_hist  = pset.getParameter<bool>("useMistagHist");

  // initialize the mistag calculation using CTEQ6L pdfs
  if (beam_energy() > 0 && _calculate_mistag) {
    char pdfname[] = "cteq6l.LHpdf";
    _mistag_calc = new MistagCalc(2*beam_energy(), pdfname);
  }

  static const double twopi = 2*3.14159265358979323846;

  _a_lim[0] = -1;
  _a_lim[1] = -_max_rap;
  _a_lim[2] = 2e-6;
  _a_lim[3] = 0;
  _a_lim[4] = _fit_win[0];
  _a_lim[5] = 0;
  _b_lim[0] = 1;
  _b_lim[1] = _max_rap;
  _b_lim[2] = _max_pt;
  _b_lim[3] = twopi;
  _b_lim[4] = _fit_win[1];
  _b_lim[5] = twopi;
  _limit_low[0] = -1;
  _limit_low[1] = -999;
  _limit_low[2] = 2e-6;
  _limit_low[3] = 0;
  _limit_low[4] = 0.1;
  _limit_low[5] = 0;
  _limit_upp[0] = 1;
  _limit_upp[1] = 999;
  _limit_upp[2] = 10*_max_pt;
  _limit_upp[3] = twopi;
  _limit_upp[4] = 1e20;
  _limit_upp[5] = twopi;
}

MistagCalc* AsymFitManager::mistag_calc() const { 
  if (_mistag_calc)
    return _mistag_calc;
  else
    throw cms::Exception("MistagCalcNotInitialized");
}

void AsymFitManager::loadParametrization(const char* cache_fn) {
  TFile* paramFile = new TFile(cache_fn, "read");
  if (!paramFile || !paramFile->IsOpen())
    throw cms::Exception("CannotOpenParamsFile") << "filename: " << cache_fn << "\n";
  TDirectory* d = (TDirectory*)paramFile->Get("AsymmetryParametrizer");

  TArrayD* arr;
  d->GetObject("mistag_pars", arr);
  for (int i = 0; i < 6; i++) mistag_pars[i] = arr->At(i);
  d->GetObject("mass_pars", arr);
  for (int i = 0; i < 7; i++) mass_pars[i] = arr->At(i);
  d->GetObject("rap_pars", arr);
  for (int i = 0; i < 5; i++) rap_pars[i] = arr->At(i);
  d->GetObject("pt_pars", arr);
  for (int i = 0; i < 5; i++) pt_pars[i] = arr->At(i);
  d->GetObject("phi_cs_pars", arr);
  for (int i = 0; i < 5; i++) phi_cs_pars[i] = arr->At(i);

  // After doing SetDirectory(0) below, we own the pointers to the
  // histograms and must delete them when done.

  h2_mistagProb = (TH2F*)d->Get("h2_mistagProb");
  h2_mistagProb->SetDirectory(0);
  h_rap_mistag_prob = (TH1F*)d->Get("h_rap_mistag_prob");
  h_rap_mistag_prob->SetDirectory(0);
  h_cos_true_mistag_prob = (TH1F*)d->Get("h_cos_true_mistag_prob");
  h_cos_true_mistag_prob->SetDirectory(0);
  h_pL_mistag_prob = (TH1F*)d->Get("h_pL_mistag_prob");
  h_pL_mistag_prob->SetDirectory(0);
  h2_pL_mistag_prob = (TH2F*)d->Get("h2_pL_mistag_prob");
  h2_pL_mistag_prob->SetDirectory(0);
  h2_cos_cs_vs_true = (TH2F*)d->Get("h2_cos_cs_vs_true");
  h2_cos_cs_vs_true->SetDirectory(0);
  h2_pTrap_mistag_prob = (TH2F*)d->Get("h2_pTrap_mistag_prob");
  h2_pTrap_mistag_prob->SetDirectory(0);

  paramFile->Close();
  delete paramFile;

  _param_cache_file = cache_fn;
}

std::ostream& operator<<(std::ostream& out, const AsymFitManager& afd) {
  out << "Mass distribution type: " << afd.mass_type() << std::endl
      << "Max rapidity: " << afd.max_rap()
      << " Max pT: " << afd.max_pt() << std::endl
      << "Fit window: " << afd.fit_win(0) << "-" << afd.fit_win(1) << std::endl
      << "Rec sigma: ";
  for (int i = 0; i < 6; i++) out << afd.rec_sigma(i) << " ";
  out << "\n";

  if (!afd.correct_mistags())
    out << "NOT including omega(...) mistag correction in pdf.\n";
  else {
    if (afd.use_mistag_hist())
      out << "Using mistag probability histogram (2D in cos_cs and rap) for mistag correction.\n";
    else if (afd.calculate_mistag())
      out << "Calculating omega(y, M) from CTEQ6L pdfs event by event for mistag correction (with omega(cos_cs) taken still from parameterization).\n";
    else
      out << "Using parameterization of omega(y, cos_cs) for mistag correction.\n";
  }

  out << "Parameterizations from file " << afd.param_cache_file() << "\n";

  out << std::setw(15) << "mistag_pars = ";
  for (int i = 0; i < 6; i++)
    out << std::setw(9) << std::setprecision(3) << afd.mistag_pars[i] << " "; 
  out << "\n" << std::setw(15) << "rap_pars = ";
  for (int i = 0; i < 5; i++)
    out << std::setw(9) << std::setprecision(3) << afd.rap_pars[i] << " "; 
  out << "\n" << std::setw(15) << "pt_pars = ";
  for (int i = 0; i < 5; i++)
    out << std::setw(9) << std::setprecision(3) << afd.pt_pars[i] << " "; 
  out << "\n" << std::setw(15) << "phi_cs_pars = ";
  for (int i = 0; i < 5; i++)
    out << std::setw(9) << std::setprecision(3) << afd.phi_cs_pars[i] << " "; 
  out << "\n" << std::setw(15) << "mass_pars = ";
  for (int i = 0; i < 7; i++)
    out << std::setw(9) << std::setprecision(3) << afd.mass_pars[i] << " "; 

  return out;
}
