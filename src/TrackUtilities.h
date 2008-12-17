#ifndef TRACKUTILITIES_H
#define TRACKUTILITIES_H

#include "DataFormats/TrackReco/interface/Track.h"

// Get the error on the inverse quantity by simple error propagation.
inline double invError(double val, double err) {
  return 1/val/val*err;
} 

// Track error wrapper functions.

inline double pError(const reco::Track* tk) {
  return invError(1/tk->p(), tk->qoverpError());
}

inline double ptError(const reco::Track* tk) {
  return tk->ptError();
}

inline double invPtError(const reco::Track* tk) {
  return invError(tk->pt(), tk->ptError());
}
  
inline double invPError(const reco::Track* tk) {
  return tk->qoverpError();
}

#endif
