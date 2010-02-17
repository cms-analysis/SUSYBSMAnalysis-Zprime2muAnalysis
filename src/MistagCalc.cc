#include <cmath>
#include <cstring>
#include <iostream>

#include "MistagCalc.h"

// the charges squared of u- and d-type quarks
const double MistagCalc::chgsq[2] = {4./9, 1./9};

MistagCalc::MistagCalc(const double sqrt_s_, char* pdfName, int subset)
  : sqrt_s(sqrt_s_) {
  lhapdf::initpdfsetbyname_(pdfName, std::strlen(pdfName));
  lhapdf::initpdf_(subset);
}

void MistagCalc::xfx(double x, double Q, double* newf) {
  lhapdf::evolvepdf_(x, Q, newf);
}

double MistagCalc::omega(const double y, const double m) {
  // Calculate the mistag probability for rapidity y and diquark
  // invariant mass m. Returns 0 on failure (i.e. unphysical
  // combination of y and m).

  // The kinematics, assuming pT = 0:
  // Parton momenta:
  //   p1 = Ebeam(x1, 0, 0,  x1)
  //   p2 = Ebeam(x1, 0, 0, -x2)
  // which implies
  //   m^2 = (p1 + p2)^2 = x1 x2 s = 4 x1 x2 Ebeam^2,
  //   y = 0.5 ln[(Ez + Pz)/(Ez - Pz)] = 0.5 ln(x1/x2);
  // solve to get x1,2 = m/sqrt(s) exp(\pm abs(y)).
  // (Only looking at abs(y) so that x1 > x2.)
  // From assuming that p_q \parallel p_Z, we mistag if the quark ends
  // up with x2 and the anti-quark with x1.

  double omega = 0;
  if (cache_.has(y,m,omega))
    return omega;

  if (m < 1e-6 || m > sqrt_s) {
    std::cerr << "Unphysical m = " << m << " in MistagCalc::omega(y,m)!" << std::endl;
    return 0;
  }
  
  const double absy = fabs(y);
  double x1, x2, f[2][13];
  x1 = m/sqrt_s*exp(absy);
  if (x1 > 1) return 0;
  x2 = m/sqrt_s*exp(-absy); // if we get to here, x2 is guaranteed < 1

  // evolve the pdfs
  xfx(x1, m, f[0]);
  xfx(x2, m, f[1]);

  double prob = 0, mistag = 0;
  for (unsigned q = 1; q <= 6; q++) {
    prob += chgsq[q%2]*(f[0][q+6]*f[1][6-q] + f[1][q+6]*f[0][6-q]);
    // from the above, x1 > x2; we mistag if the quark has x2 and the
    // anti-quark has x1.
    mistag += chgsq[q%2]*f[1][q+6]*f[0][6-q];
  }

  omega = prob == 0 ? 0 : mistag/prob;
  cache_.put(y,m,omega);
  return omega;
}
