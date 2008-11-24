#ifndef DILEPTONUTILITIES_H
#define DILEPTONUTILITIES_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "FWCore/Framework/interface/Event.h"

// Dilepton utility functions

// Sort a dilepton collection by decreasing invariant mass, and then
// prune the collection: if more than one dilepton was formed,
// accept only those containing distinct leptons. If a lepton is
// shared by a pair of dileptons, keep the dilepton with the higher
// invariant mass.
void removeDileptonOverlap(reco::CompositeCandidateCollection& dileptons,
			   const bool debug=false);

// Take the input dileptons in dils (e.g. the output of CandCombiner)
// which have only the combinatorics done, and put only the dileptons
// that correspond to generator-level resonances into newDils.
void genDileptonsOnly(const reco::CompositeCandidateCollection& dils,
		      reco::CompositeCandidateCollection& newDils,
		      const bool debug=false);

// Take the input dileptons in dils (e.g. the output of CandCombiner)
// which have only the combinatorics done, and apply the
// analysis-level cuts (specified by a bitmask) to them. If both PDG
// ids are nonzero, cut out the dileptons that are not made up of the
// requested leptons. (Useful for separating mu+mu+ from mu-mu- in the
// output of CandCombiner.)
void cutDileptons(const edm::Event& event,
		  const reco::CompositeCandidateCollection& dils,
		  reco::CompositeCandidateCollection& newDils,
		  unsigned cuts, int pdgId1=0, int pdgId2=0,
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

#endif
