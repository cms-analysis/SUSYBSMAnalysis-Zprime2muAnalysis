#ifndef GENERATORUTILITIES_H
#define GENERATORUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Particle.h"

// Generator-level utility functions.

// Work around Cartesian v. polar coordinates for LorentzVectors for now.
void SetP4M(reco::Particle::LorentzVector& v,
	    double pt, double phi, double p,
	    double theta, double mass);

// Return the lepton cand's mother. If cand's mother is the same
// lepton but before brem, go up the decay chain (what we call the
// "non-brem mother"). This is mainly to avoid stopping at leptons
// in documentation lines, which are declared to be ancestors of
// muons produced in hard interaction (i.e., ancestors of
// themselves).
const reco::Candidate* mother(const reco::CandidateBaseRef& cand);

// Return the PDG id of the non-brem mother of cand. If the mother
// pointer isn't valid, return pdgId = 0.
int motherId(const reco::CandidateBaseRef& cand);
  
// Return the PDG id of the grandmother of cand (i.e. a quark or
// antiquark, or a gluon). If the mother or grandmother pointer isn't
// valid, return pdgId = 0.
int grandmotherId(const reco::CandidateBaseRef& cand);

// If cand1 and cand2 have the same non-brem mothers, return a
// pointer to the mother candidate, else return null.
const reco::Candidate* sameMother(const reco::CandidateBaseRef& cand1,
				  const reco::CandidateBaseRef& cand2);

#endif
