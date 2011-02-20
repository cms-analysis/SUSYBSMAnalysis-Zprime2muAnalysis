#ifndef AsymFitManager_h
#define AsymFitManager_h

#include <ostream>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"

namespace edm {
  class ParameterSet;
}

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/MistagCalc.h"

// AsymFitManager is designed to protect some constants that are used
// by the fitting routines that used to be compiled global constants,
// but now are read in from a cfg file. With accessors, this should
// protect them.  We would put all the fitting functions in a class
// together with these constants, but eh...

// AsymFitManager also is responsible for switching the cos_cs pdf,
// the asym_3_PDF function pointer in AsymFunctions.C; since we need
// to be able to fit to different forms (for true asymmetry fits, for
// gravitons, etc).

// Also, we can still forget to set all the values before calling the
// accessors, but this is at least a little safer, code-wise (we won't
// later be able to easily overwrite the values).

class AsymFitManager {
 public:
  AsymFitManager() 
    : _mistag_calc(0),
      _param_cache_file("none yet!"),
      h2_mistagProb(0),
      h_rap_mistag_prob(0),
      h_cos_true_mistag_prob(0),
      h_pL_mistag_prob(0),
      h2_pL_mistag_prob(0),
      h2_cos_cs_vs_true(0),
      h2_pTrap_mistag_prob(0)
  {}

  ~AsymFitManager() {
    delete _mistag_calc;
    delete h2_mistagProb;
    delete h_rap_mistag_prob;
    delete h_cos_true_mistag_prob;
    delete h_pL_mistag_prob;
    delete h2_pL_mistag_prob;
    delete h2_cos_cs_vs_true;
    delete h2_pTrap_mistag_prob;
  }

  void setConstants(const edm::ParameterSet& pset);
  void loadParametrization(const char* cache_fn);

  // A small number to return instead of zero in pdf calculation.
  double epsilon() const { return 1.e-6; }
  // Roughly model the acceptance with a cut in eta at 2.4.
  double eta_lim() const { return 2.4; }
  // Can cut separately on negative and positive muons.
  double mup_eta_lim_lo() const { return -2.4; }
  double mup_eta_lim_hi() const { return  2.4; }
  double mum_eta_lim_lo() const { return -2.4; }
  double mum_eta_lim_hi() const { return  2.4; }
  double pt_min() const { return 0.; } // used 20 for gravitons
  double mup_pt_min() const { return 0.; }
  double mum_pt_min() const { return 0.; }
  bool debug() const { return false; }

  bool doing_electrons() const { return _doing_electrons; }
  double lepton_mass() const { return _lepton_mass; }
  int mass_type() const { return _mass_type; } 
  double max_rap() const { return _max_rap; }
  double max_pt() const { return _max_pt; } 
  double fit_win(unsigned int i) const { return _fit_win[i]; }
  double rec_sigma(unsigned int i) const { return _rec_sigma[i]; }
  double a_lim(unsigned int i) const { return _a_lim[i]; }
  double b_lim(unsigned int i) const { return _b_lim[i]; }
  const double* a_lim_arr() const { return _a_lim; }
  const double* b_lim_arr() const { return _b_lim; }
  double limit_low(unsigned int i) const { return _limit_low[i]; }
  double limit_upp(unsigned int i) const { return _limit_upp[i]; }
  double peak_mass() const { return _peak_mass; }
  double beam_energy() const { return 3500.; }

  MistagCalc* mistag_calc() const;

  enum RECSIGMA { SIGMA_COSCS, SIGMA_RAP, SIGMA_PT,
		  SIGMA_PHI, SIGMA_MASS, SIGMA_PHICS };

  double recSigmaCosCS() const { return _rec_sigma[SIGMA_COSCS]; }
  double recSigmaRap() const { return _rec_sigma[SIGMA_RAP]; }
  double recSigmaPt() const { return _rec_sigma[SIGMA_PT]; }
  double recSigmaPhi() const { return _rec_sigma[SIGMA_PHI]; }
  double recSigmaMass() const { return _rec_sigma[SIGMA_MASS]; }
  double recSigmaPhiCS() const { return _rec_sigma[SIGMA_PHICS]; }

  bool correct_mistags()  const { return _correct_mistags;  }
  bool calculate_mistag() const { return _calculate_mistag; }
  bool use_mistag_hist()  const { return _use_mistag_hist;  }

  const std::string& param_cache_file() const { return _param_cache_file; }

  enum MASS_TYPE { MASS_EXP=1, MASS_LOR, MASS_LOREXP };

 private:
  MistagCalc* _mistag_calc;

  bool _doing_electrons;
  double _lepton_mass;
  // mass_type values: 1 = falling exponential, 2 = lorentzian peak, 3 = 1+2
  int _mass_type;
  double _max_rap;
  double _max_pt;
  std::vector<double> _fit_win;
  double _peak_mass;
  // These sigma are the reconstructed - generated for dilepton pt, rap, mass, 
  // phi, cos_cs, and phi_cs.  They are taken from DY generated above 1 TeV.
  // order: cos_cs, rapidity, pT, phi, mass, phi_cs (see above enum)
  std::vector<double> _rec_sigma;
  // These are the limits of integration used for normalized the smeared 
  // PDF.  Integration must be done over all possible reconstructed values
  // (not std::vectors so we can pass to a function that doesn't expect one)
  double _a_lim[6];
  double _b_lim[6];
  double _limit_low[6];
  double _limit_upp[6];

  bool _correct_mistags;
  bool _calculate_mistag;
  bool _use_mistag_hist;

  std::string _param_cache_file;

public:
  double mistag_pars[6];
  double mass_pars[7];
  double rap_pars[5];
  double pt_pars[5];
  double phi_cs_pars[5];
  
  TH2F* h2_mistagProb;
  TH1F* h_rap_mistag_prob;
  TH1F* h_cos_true_mistag_prob;
  TH1F* h_pL_mistag_prob;
  TH2F* h2_pL_mistag_prob;
  TH2F* h2_cos_cs_vs_true;
  TH2F* h2_pTrap_mistag_prob;
};

std::ostream& operator<<(std::ostream& out, const AsymFitManager& afd);

#endif
