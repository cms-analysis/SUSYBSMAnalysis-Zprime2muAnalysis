#ifndef MISTAGCALC_H
#define MISTAGCALC_H

#include "MistagCache.h"

namespace lhapdf {
  extern "C" {
    void initpdfsetbyname_ (char *, int len);
    void initpdf_(int &);
    void evolvepdf_(double &, double &, double *);
  }
}

class MistagCalc {
public:
  MistagCalc(const double pb, char* pdfName, int subset=0);
  ~MistagCalc();

  static const double chgsq[2];

  void xfx(double x, double Q2, double* newf);
  double omega(const double y, const double m);

private:
  // the center-of-mass energy in GeV
  double sqrt_s;
  // a cache to store calculated values in
  mtagcache* cache;
};

#endif
