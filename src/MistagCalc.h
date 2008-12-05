#ifndef MISTAGCALC_H
#define MISTAGCALC_H

#include <map>

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

  static const double chgsq[2];

  void xfx(double x, double Q2, double* newf);
  double omega(const double y, const double m);

private:
  class cache {
  public:
    bool has(double y, double m, double& omega) {
      cachetype::const_iterator it = cache_.find(std::make_pair(y,m)); 
      if (it != cache_.end()) {
	omega = it->second;
	return true;
      }
      return false;
    }
  
    void put(double y, double m, const double omega) {
      cache_[std::make_pair(y,m)] = omega;
    }
  
  private:
    typedef std::map<std::pair<double, double>, double> cachetype;
    cachetype cache_;
  };

  // the center-of-mass energy in GeV
  double sqrt_s;
  // a cache to store calculated values in
  cache cache_;
};

#endif
