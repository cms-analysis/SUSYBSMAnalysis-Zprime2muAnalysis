#ifndef DILEPTONUTILITIES_H
#define DILEPTONUTILITIES_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

// Sorting functor for dileptons. Using std::sort results in dileptons
// reverse-sorted by invariant mass.
struct reverse_mass_sort {
  bool operator()(const reco::CompositeCandidate& lhs, const reco::CompositeCandidate& rhs) {
    return lhs.mass() > rhs.mass();
  }
};

// Sort a dilepton collection by decreasing invariant mass, and then
// prune the collection: if more than one dilepton was formed,
// accept only those containing distinct leptons. If a lepton is
// shared by a pair of dileptons, keep the dilepton with the higher
// invariant mass.
void removeDileptonOverlap(reco::CompositeCandidateCollection& dileptons,
			   const bool debug=false);

// Count the number of daughters the dilepton has in the specified
// acceptance in eta.
int numDaughtersInAcc(const reco::CompositeCandidate& dil,
		      const double etaCut=2.4);

// Return a reference to the ith daughter lepton of the dilepton, or
// an invalid reference if i is out of bounds.
const reco::CandidateBaseRef
dileptonDaughter(const reco::CompositeCandidate& dil,
		 const unsigned i);

// Return a reference to the daughter lepton of the dilepton with
// specified charge (if it is a same-sign dilepton, will return the
// first one found), or else an invalid reference if not found.
// (Code using this method should check the ref for validity, since
// for same-sign dileptons it is easy to not find a daughter lepton
// of the wrong sign.)
const reco::CandidateBaseRef
dileptonDaughterByCharge(const reco::CompositeCandidate& dil,
			 const int charge);

// Get the momentum four-vector of the composite candidate plus the
// leptons' closest photons (what we've in the past called the
// "resonance").
reco::Particle::LorentzVector resonanceP4(const pat::CompositeCandidate& cand);

#endif
