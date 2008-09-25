#include <cmath>

#include "Utilities.h"

namespace Utilities {
  std::string binary(unsigned n, bool pad) {
    std::string s;
    while (n) {
      if (n & 1) s = "1" + s;
      else       s = "0" + s;
      n >>= 1;
    }
    if (pad)
      for (unsigned i = 0; i < 8 - s.size()%8; ++i)
	s = "0" + s;
    return s;
  }

  double delta_phi(double phi1, double phi2) {
    static const double pi = 3.14159265;
    double d = phi1 - phi2;
    if      (d >  pi) d -= 2*pi;
    else if (d < -pi) d += 2*pi;
    return d;
  }

  double delta_r(double eta1, double eta2, double phi1, double phi2) {
    double dp = delta_phi(phi1, phi2);
    double de = eta1 - eta2;
    return sqrt(dp*dp + de*de);
  }

  bool HistogramIncreasingIntegral(const histo_bkg& a, const histo_bkg& b) {
    return a.first->Integral(a.first->FindBin(200), a.first->FindBin(1e9)) <
           b.first->Integral(b.first->FindBin(200), b.first->FindBin(1e9));
  }

  bool HistogramDecreasingMaximum(const histo_bkg& a, const histo_bkg& b) {
    return a.first->GetMaximum() > b.first->GetMaximum();
  }
}
