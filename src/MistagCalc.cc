#include <cmath>
#include <iostream>

#include "MistagCalc.h"
#include "MistagCache.h"

bool operator<(const mtagdata& l, const mtagdata& r) {
  if (l.y < r.y) return true;
  else if (l.y > r.y) return false;
  else {
    if (l.m < r.m) return true;
    else if (l.m > r.m) return false;
    else {
      // they are equal
      return false;
    }
  }
}

// the charges squared of u- and d-type quarks
const double MistagCalc::chgsq[2] = {4./9, 1./9};

MistagCalc::MistagCalc(const double sqrt_s_, char* pdfName, int subset)
  : sqrt_s(sqrt_s_) {
  lhapdf::initpdfsetbyname_(pdfName, strlen(pdfName));
  lhapdf::initpdf_(subset);
  cache = new mtagcache;
}

MistagCalc::~MistagCalc() {
  delete cache;
}

void MistagCalc::xfx(double x, double Q, double* newf) {
  lhapdf::evolvepdf_(x, Q, newf);
}

double MistagCalc::omega(const double y, const double m) {
  // Calculate the mistag probability for rapidity y and diquark
  // invariant mass m. Returns 0 on failure (i.e. unphysical
  // combination of y and m).

  double omega = 0;
  if (cache->has(y,m,omega))
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
  cache->put(y,m,omega);
  return omega;
}
