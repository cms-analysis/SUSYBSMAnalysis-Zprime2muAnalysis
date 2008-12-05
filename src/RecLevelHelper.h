#ifndef RECLEVELHELPER_H
#define RECLEVELHELPER_H

#include <map>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// details about the number of rec levels stored, their names, etc.
const int TRIG_LEVELS = 4; // Includes lgen == 0 as a placeholder.
enum RecLevel { lgen, l1, l2, l3, lgmr, ltk, lfms, lpmr, lbest };
const int MAX_LEVELS = lbest+1;

// The names of the defined rec levels, as well as shortened versions.
const std::string levelNames[MAX_LEVELS] = {
  "Gen", " L1", " L2", " L3", "GMR", "Tracker-only", "TPFMS", "PMR", "OPT"
};
const std::string levelNamesShort[MAX_LEVELS] = {
  "GN", "L1", "L2", "L3", "GR", "TK", "FS", "PR", "OP"
};

// Perform a sanity check on the rec level passed, throwing an
// exception if it is out of range.
void checkRecLevel(const int level, const char* name);

// Return either the short or the regular version of rec level's
// name.
const std::string& levelName(const int rec, bool shortVersion=false);

class RecLevelHelper {
 public:
  // Set up per-lifetime things, including storing the collection
  // InputTags in the right order.
  void init(const edm::ParameterSet& config);

  // Set up per-event things, such as the recLevelMap, and match maps.
  void initEvent(const edm::Event& event);

  // Using getLeptonsView(), get a vector of refs to the lepton
  // candidates at the requested rec level, storing them in leps.
  bool getLeptons(const edm::Event& event, int level,
		  reco::CandidateBaseRefVector& leps);

  // Same, but for the dileptons.
  bool getDileptons(const edm::Event& event, int level,
		    reco::CompositeCandidateCollection& coll);

  // Return whether the rec level has a collection in the event.
  bool recLevelOkay(const edm::Event& event, int level);

  // Return a unique id for the lepton cand (currently implemented by
  // the reference's index into the collection).
  int id(const reco::CandidateRef& cand) const;
  int id(const reco::CandidateBaseRef& cand) const;

  // Translate the Ref's product id to one of our rec levels, using
  // the cached map. If the rec level turns out to be lbest, return
  // the original rec level using the stored map.
  int recLevel(const reco::CandidateBaseRef& cand) const;

  // Get the rec level for a dilepton, making sure that either all the
  // daughter leptons have the same rec level, or else returning lbest,
  // (since a "best" dilepton can be made up of leptons at different rec
  // levels).
  int recLevel(const reco::CompositeCandidate& cand) const;

  // Find the original rec level of cand; useful for "best" leptons to
  // find which of ltk, lfms, lpmr the cocktail chose.
  int originalRecLevel(const reco::CandidateBaseRef& cand) const;

  // Get the four-vector of the closest photon found for cand.
  reco::Particle::LorentzVector
    closestPhoton(const reco::CandidateBaseRef& cand) const;

  // Return the already-matched-in-deltaR closest MC truth
  // lepton. Returns an invalid ref if not found.
  const reco::CandidateBaseRef&
    matchGenLepton(const reco::CandidateBaseRef& lep) const;

  // Get the id of lep's closest match (in delta R) at another rec level.
  int genMatchId(const reco::CandidateBaseRef& lep) const
  { return id(matchGenLepton(lep)); }

  // Get the lepton of the lepton at the other offline level 
  const reco::CandidateBaseRef&
    matchOfflineLepton(const reco::CandidateBaseRef& lep,
		       const int level) const;


 private:
  // The arrays of InputTags which maps rec level number to the EDM
  // collections.
  edm::InputTag lepInputs[MAX_LEVELS];
  edm::InputTag dilInputs[MAX_LEVELS];

  // Flag whether we've warned about a missing collection at each rec
  // level for each event. (Set to prevent double warnings when
  // getLeptonsView() gets called twice.)
  mutable bool warned[MAX_LEVELS];

  // Map product IDs from CandidateBaseRef to whatever rec level we
  // choose.
  std::map<int,int> recLevelMap;
  
  // Map lepton CandidateBaseRefs to photon ones.
  reco::CandViewMatchMap photonMatchMap[MAX_LEVELS];

  // Store an AssociationMap for each rec level which maps leptons to
  // leptons at other rec levels, matching by closest in delta R. (We
  // are wasteful, not using entries for which the rec levels are
  // equal, but storage is easier this way.)
  reco::CandViewMatchMap genMatchMap[MAX_LEVELS];

  // An invalid CandidateBaseRef so that some of the methods above can
  // return by reference and still be able to return a null reference.
  reco::CandidateBaseRef invalidRef;

  // Store the object which maps EDM collections (Views, or members of
  // collections) to our defined rec level numbers (in a sense, the
  // inverse operation to getLeptonsView() above).
  void storeRecLevelMap(const edm::Event& event);

  // Retrieve a single match map from the event by name, and return it
  // by reference. If the match map does not exist, print a warning
  // but do not throw an exception.
  void storeMatchMap(const edm::Event& event, const std::string& mapName,
		     reco::CandViewMatchMap& map) const;

  // Store the match maps for this event: the matching between all
  // pairs of rec levels ("closest" matching), the matching between
  // all global fits ("seed" matching), and closest photon matches.
  void storeMatchMaps(const edm::Event& event);

  // Using the vector of InputTags, get a collection reference at the
  // requested rec level.
  template <typename T>
    bool get(const edm::Event& event, int recLevel,
	     const edm::InputTag& tag, T& coll) const;
};

#endif
