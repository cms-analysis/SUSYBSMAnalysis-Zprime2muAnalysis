#ifndef TRACKUTILITIES_H
#define TRACKUTILITIES_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Track utility functions

// There is no pError() in Track/TrackBase; calculate error on p
// ourselves from error on qoverp().
double pError(const reco::Track* cand);
// There is a ptError() in Track/TrackBase but for symmetry with the
// above let's have another.
double ptError(const reco::Track* cand);
  
// Propagate inverse errors...
inline double invError(double val, double err) { return 1/val/val*err; } 
// ... to 1/pT
double invPtError(const reco::Track* tk);
// ... and 1/p
double invPError(const reco::Track* tk);

// Return whether cand1 and cand2 are "close"; instead of a circle
// in eta-phi space, we look at a square .5 on a side.
bool matchTracks(const reco::CandidateBaseRef& cand1,
		 const reco::CandidateBaseRef& cand2);

#endif
