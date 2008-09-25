#ifndef Utilities_h
#define Utilities_h

#include <string>

#include "TH1F.h"

#include "Parameters.h"

namespace Utilities {
  // Return a string with n in binary representation. If desired, pad
  // the string to the nearest multiple of 8 bits.
  std::string binary(unsigned n, bool pad=true);

  // Compute delta_phi \in (-pi, pi] and delta R = sqrt(delta_phi^2 + delta_eta^2).
  double delta_phi(double phi1, double phi2);
  double delta_r(double eta1, double eta2, double phi1, double phi2);

  // Associate a histogram to one of our background id codes.
  typedef std::pair<TH1F*, Parameters::bkg_id_type> histo_bkg;
  
  // Helper function for sorting mass histograms by increasing integrals above 200 (GeV).
  bool HistogramIncreasingIntegral(const histo_bkg& a, const histo_bkg& b);

  // Helper function for sorting histograms by decreasing maxima.
  bool HistogramDecreasingMaximum(const histo_bkg& a, const histo_bkg& b);
}

#endif

