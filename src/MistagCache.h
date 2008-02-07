#ifndef MISTAGCACHE_H
#define MISTAGCACHE_H

#include <map>

struct mtagdata {
  mtagdata(double y_, double m_) : y(y_), m(m_) {}

  double y;
  double m;

  friend bool operator<(const mtagdata& l, const mtagdata& r);
};

class mtagcache {
 public:
  mtagcache() {
    hit = miss = 0;
  }
  ~mtagcache() {
    printf("mtagcache stats: hit %i miss %i total %i\n", hit, miss, hit+miss);
  }
    
  bool has(double y, double m, double& omega) {
    mtagdata thismt(y,m);
    std::map<mtagdata, double>::const_iterator it = cache.find(thismt); 
    if (it != cache.end()) {
      omega = it->second;
      hit++;
      return true;
    }
    miss++;
    return false;
  }
  
  void put(double y, double m, const double omega) {
    mtagdata thismt(y,m);
    cache[thismt] = omega;
  }
  
 private:  
  std::map <mtagdata, double> cache;
  int hit;
  int miss;
};

#endif
