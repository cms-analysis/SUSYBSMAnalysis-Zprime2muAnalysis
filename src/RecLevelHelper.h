#ifndef RECLEVELHELPER_H
#define RECLEVELHELPER_H

#include <map>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// details about the number of rec levels stored, their names, etc.
const int NUM_REC_LEVELS = 4;
const int MAX_LEVELS = 8;
enum RecLevel { lgen, l1, l2, l3, lgmr, ltk, lfms, lpmr, lbest };
const std::string str_level[MAX_LEVELS+1] = {
  "Gen", " L1", " L2", " L3", "GMR", "Tracker-only", "TPFMS", "PMR", "OPT"
};
const std::string str_level_short[MAX_LEVELS+1] = {
  "GN", "L1", "L2", "L3", "GR", "TK", "FS", "PR", "BS"
};

class RecLevelHelper {
 public:
  bool isDoingElectrons() { return doingElectrons; }
  
  // Set up per-lifetime things, including storing the collection
  // InputTags in the right order (padding if doingElectrons).
  void init(const edm::ParameterSet& config);

  // Set up per-event things, such as the recLevelMap, and match maps.
  void initEvent(const edm::Event& event);

  // Using the vector of InputTags, get a view to the candidates at
  // the requested rec level.
  bool getView(const edm::Event& event,
	       int recLevel, edm::View<reco::Candidate>& view);

  // Helper enums and function to build the match map name for storing
  // in the event. E.g. storing a closest match map from L3 to GMR
  // will get the branch name *_leptonMatches_closestL3GR_*.
  enum MatchType { PHOTON, CLOSEST, SEED };
  std::string makeMatchMapName(MatchType mtype, const int irec,
			       const int jrec=-1) const;

  // Translate the Ref's product id to one of our rec levels, using the
  // cached map.
  int recLevel(const reco::CandidateBaseRef& cand) const;

  // Perform a sanity check on the rec level passed, throwing an
  // exception if it is out of range.
  void checkRecLevel(const int level, const char* name) const;

  // Get the four-vector of the closest photon found for cand.
  reco::Particle::LorentzVector
    closestPhoton(const reco::CandidateBaseRef& cand) const;

  // Get the seed index (i.e. the index into the stand-alone muon
  // collection) of the candidate.
  int seedIndex(const reco::CandidateBaseRef& cand) const;

  // Search for the lepton at the specified rec level which is either
  // the closest or same-seed match to the lepton specified, returning
  // an invalid reference if not found (if the id was = -999 or an
  // invalid rec level). whichMatch = 0 for closest match, 1 for
  // same-seed match, and -1 to pick the "best" match, i.e. same-seed if
  // available, and closest if not.  Not meant to be called directly;
  // use one of the three methods closestLepton, sameSeedLepton, and
  // matchedLepton.
  const reco::CandidateBaseRef&
    matchLepton(const reco::CandidateBaseRef& lep,
		const int level,
		int whichMatch) const;

  // Const access to lep's closest match (in delta R) at another rec level.
  const reco::CandidateBaseRef&
    closestLepton(const reco::CandidateBaseRef& lep,
		  const int level) const;

  // Const access to lep's same-seed match at another rec level.
  const reco::CandidateBaseRef&
    sameSeedLepton(const reco::CandidateBaseRef& lep,
		   const int level) const;

  // Const access to lep's "best" match (defined above) at another rec
  // level.
  const reco::CandidateBaseRef&
    matchedLepton(const reco::CandidateBaseRef& lep,
		  const int level) const;

 private:
  // Flag stating whether we're running in electron mode, which 
  // ignores all collections except gen and gmr.
  bool doingElectrons;

  // The vector of InputTags which maps rec level number to the EDM
  // collection.
  std::vector<edm::InputTag> inputs;

  // Flag whether we've warned about a missing collection at each rec
  // level for each event. (Set to prevent double warnings when
  // getView() gets called twice.)
  bool warned[MAX_LEVELS];

  // Map product IDs from CandidateBaseRef to whatever rec level we
  // choose.
  std::map<int,int> recLevelMap;
  
  // Map lepton CandidateBaseRefs to photon ones.
  reco::CandViewMatchMap photonMatchMap[MAX_LEVELS];

  // Store seed indices for global fits (indexed in the vector by
  // id(cand)).
  std::vector<int> seedIndices[MAX_LEVELS];

  // Store an AssociationMap for each rec level which maps leptons to
  // leptons at other rec levels, matching by closest in delta R. (We
  // are wasteful, not using entries for which the rec levels are
  // equal, but storage is easier this way.)
  reco::CandViewMatchMap closestMatchMap[MAX_LEVELS][MAX_LEVELS];

  // Store an AssociationMap for each rec level of global fits
  // (i.e. L3 and above), matching by "seed" lepton. (Currently only
  // meaningful for muons, and we are even more wasteful than before,
  // in addition not using entries for which rec levels are less than
  // l3.)
  reco::CandViewMatchMap seedMatchMap[MAX_LEVELS][MAX_LEVELS];

  // An invalid CandidateBaseRef so that some of the methods above can
  // return by reference and still be able to return a null reference.
  reco::CandidateBaseRef invalidRef;

  // Store the object which maps EDM collections (Views, or members of
  // collections) to our defined rec level numbers (in a sense, the
  // inverse operation to getView() above).
  void storeRecLevelMap(const edm::Event& event);

  // Store the match maps for this event: the matching between all
  // pairs of rec levels ("closest" matching), the matching between
  // all global fits ("seed" matching), and closest photon matches.
  void storeMatchMaps(const edm::Event& event);
};

#endif
