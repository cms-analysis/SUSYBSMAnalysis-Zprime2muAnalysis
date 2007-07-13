#ifndef ASYMFITMANAGER_H
#define ASYMFITMANAGER_H

#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

// forward declare for AsymFunctions.h
class AsymFitManager;

#include "AsymFunctions.h"

// AsymFitManager is designed to protect some constants that are used
// by the fitting routines that used to be compiled global constants,
// but now are read in from a cfg file. With accessors, this should
// protect them.  We would put all the fitting functions in a class
// together with these constants, but ROOT cannot use a member
// function with its TF1 class for fitting!

// AsymFitManager also is responsible for switching the cos_cs pdf,
// the asym_3_PDF function pointer in AsymFunctions.C; since we need
// to be able to fit to different forms (for true asymmetry fits, for
// gravitons, etc).

// Issue: will the compiler inline the accessors correctly, so that
// the calls are quick? The accessor will be called inside the
// high-level fit functions, e.g. recAsym2D, so they shouldn't be
// called *too* many times.  Explicitly inline them, just in
// case, but does this help?

// Also, we can still forget to set all the values before calling the
// accessors, but this is at least a little safer, code-wise (we won't
// later be able to easily overwrite the values).

// JMTBAD ROOT 5.16 now can use a class member function! What version
// of CMSSW uses this ROOT version?

class AsymFitManager {
 public:
  void SetConstants(const edm::ParameterSet& pset) {
    _mass_type = pset.getParameter<int>("massDistType");
    _max_rap = pset.getParameter<double>("maxRapidity");
    _max_pt = pset.getParameter<double>("maxPt");
    _fit_win = pset.getParameter<std::vector<double> >("fitWindow");
    _gen_win = pset.getParameter<std::vector<double> >("genWindow");
    _rec_sigma = pset.getParameter<std::vector<double> >("recSigma");

    // JMTBAD *sigh*
    _a_lim[0] = -1;
    _a_lim[1] = -_max_rap;
    _a_lim[2] = 2e-6;
    _a_lim[3] = 0;
    _a_lim[4] = _fit_win[0];
    _a_lim[5] = 0;
    _b_lim[0] = 1;
    _b_lim[1] = _max_rap;
    _b_lim[2] = _max_pt;
    _b_lim[3] = 2*TMath::Pi();
    _b_lim[4] = _fit_win[1];
    _b_lim[5] = 2*TMath::Pi();
    _limit_low[0] = -1;
    _limit_low[1] = -999;
    _limit_low[2] = 2e-6;
    _limit_low[3] = 0;
    _limit_low[4] = 0.1;
    _limit_low[5] = 0;
    _limit_upp[0] = 1;
    _limit_upp[1] = 999;
    _limit_upp[2] = 10*_max_pt;
    _limit_upp[3] = 2*TMath::Pi();
    _limit_upp[4] = 1e20;
    _limit_upp[5] = 2*TMath::Pi();

    // the default is already initialized in AsymFunctions.C,
    // but do it here to get everything synchronized
    setPDF(ASYM);
  }

  // accessors for the constants (no bounds checking!)
  int mass_type() const { return _mass_type; } 
  double max_rap() const { return _max_rap; }
  double max_pt() const { return _max_pt; } 
  double fit_win(unsigned int i) const { return _fit_win[i]; }
  double gen_win(unsigned int i) const { return _gen_win[i]; }
  double rec_sigma(unsigned int i) const { return _rec_sigma[i]; }
  double a_lim(unsigned int i) const { return _a_lim[i]; }
  double b_lim(unsigned int i) const { return _b_lim[i]; }
  const double* a_lim_arr() const { return _a_lim; }
  const double* b_lim_arr() const { return _b_lim; }
  double limit_low(unsigned int i) const { return _limit_low[i]; }
  double limit_upp(unsigned int i) const { return _limit_upp[i]; }

  // the currently available pdfs
  enum PDFTYPE { ASYM, GRAV, GRAVTH };

  // update the function pointer to the one specified by _type
  // legitimate values are defined in an enum above
  // throws if _type is unrecognized
  void setPDF(PDFTYPE _type) {
    //    extern double (*asym_3_PDF)(double *x, double *par);
    pdfType = _type;
    switch (pdfType) {
    case ASYM: asym_3_PDF = asym_3_PDF_real; break;
    case GRAV: asym_3_PDF = GravitonCos_3_PDF; break;
    case GRAVTH: asym_3_PDF = GravitonCos_th_PDF; break;
    default: 
      throw cms::Exception("AsymFitManager") 
	<< "unrecognized PDF type " << _type;
    }
  }

  // return the current pdf's name; useful in pretty printing
  const char* getPDFName() const {
    static const char* PDFNAMES[3] = {
      "asym_3_PDF",
      "GravitonCos_3_PDF",
      "GravitonCos_th_PDF"
    };
    return PDFNAMES[pdfType];
  }

 private:
  // mass_type values: 1 = falling exponential, 2 = lorentzian peak, 3 = 1+2
  int _mass_type;
  double _max_rap;
  double _max_pt;
  std::vector<double> _fit_win;
  std::vector<double> _gen_win;
  // These sigma are the reconstructed - generated for dilepton pt, rap, mass, 
  // phi, cos_cs, and phi_cs.  They are taken from DY generated above 1 TeV.
  std::vector<double> _rec_sigma;
  // These are the limits of integration used for normalized the smeared 
  // PDF.  Integration must be done over all possible reconstructed values
  // (not std::vectors so we can pass to a function that doesn't expect one)
  double _a_lim[6];
  double _b_lim[6];
  double _limit_low[6];
  double _limit_upp[6];

  PDFTYPE pdfType;
};

std::ostream& operator<<(std::ostream& out, const AsymFitManager& afd);

#endif // ASYMFITMANAGER_H
