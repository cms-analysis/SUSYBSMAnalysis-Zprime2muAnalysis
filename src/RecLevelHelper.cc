#include <string>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

using namespace std;

void checkRecLevel(const int level, const char* name) {
  if (level < 0 || level > MAX_LEVELS)
    throw cms::Exception(name)
      << "invalid level " << level << " in " << name << "!\n";
}

bool isCocktailLevel(const int level) {
  checkRecLevel(level, "isCocktail");
  return level >= lOP && level < MAX_LEVELS;
}

const std::string& levelName(const int rec, bool shortVersion) {
  checkRecLevel(rec, "levelName");
  if (shortVersion) return levelNamesShort[rec];
  else return levelNames[rec];
}

void RecLevelHelper::init(const edm::ParameterSet& config) {
  bestRecLevel = config.getParameter<int>("bestRecLevel");

  for (int i = 0; i < MAX_LEVELS; i++) {
    string tagname = "leptons" + levelName(i);
    lepInputs[i] = config.exists(tagname) ? config.getParameter<edm::InputTag>(tagname) : edm::InputTag();
    tagname = "di" + tagname;
    dilInputs[i] = config.exists(tagname) ? config.getParameter<edm::InputTag>(tagname) : edm::InputTag();
  }
}

void RecLevelHelper::initEvent(const edm::Event& event) {
  for (int rec = 0; rec < MAX_LEVELS; rec++)
    warned[rec] = false;
  storeRecLevelMap(event);
  storeMatchMaps(event);
}

template <typename T>
bool RecLevelHelper::get(const edm::Event& event, int level,
			 const edm::InputTag& tag, T& coll) const {
  edm::Handle<T> handle;
  event.getByLabel(tag, handle);

  if (handle.failedToGet()) {
    if (!warned[level]) {
      string inp = tag.encode();
      // Don't bother to warn about collections that are supposed to be missing.
      if (inp != ":" && inp != "") 
	edm::LogWarning("initEvent")
	  << "No event collection " << inp
	  << " found at rec level " << level << "; skipping";
      warned[level] = true;
    }

    return false;
  }

  coll = *handle;
  return true;
}

bool RecLevelHelper::getLeptons(const edm::Event& event, int level,
				reco::CandidateBaseRefVector& leps) {
  edm::View<reco::Candidate> view;
  if (!get(event, level, lepInputs[level], view))
    return false;

  // Store refs to the leptons.
  leps.clear();
  for (unsigned ilep = 0; ilep < view.size(); ++ilep)
    leps.push_back(view.refAt(ilep));

  return true;
}

bool RecLevelHelper::getDileptons(const edm::Event& event, int level,
				  reco::CompositeCandidateCollection& dils) {
  return get(event, level, dilInputs[level], dils);
}  

bool RecLevelHelper::recLevelOkay(const edm::Event& event, int level) {
  edm::View<reco::Candidate> view;
  return get(event, level, lepInputs[level], view);
}
  
void RecLevelHelper::storeRecLevelMap(const edm::Event& event) {
  // We shouldn't have to re-store this every event...
  recLevelMap.clear();
  for (int rec = 0; rec < MAX_LEVELS; rec++) {
    edm::View<reco::Candidate> view;
    get(event, rec, lepInputs[rec], view);

    // Cache the product ids for each collection so we can look up the
    // rec level from the CandidateBaseRef.
    int key = int(view.id().id());
    if (recLevelMap.find(key) == recLevelMap.end())
      recLevelMap[key] = rec;

    // So that we can figure out from which rec level a cocktail muon
    // originated, store a map of the product ids of the Track objects
    // for each of the possible constituent muons (e.g. rec level lGR
    // should map to the product id for the recoTracks_globalMuons_*_*
    // branch, etc.) We can store this in the same map as the above
    // since the product ids are guaranteed to be unique for each
    // branch (at least on a per-event basis).
    if (!isCocktailLevel(rec)) {
      // Get the track from the first muon in the collection for this
      // rec level. If there are no muons for this rec level, then we
      // won't be matching to this rec level this event anyway.
      if (view.size() > 0) {
	const reco::Muon* mu = toConcretePtr<reco::Muon>(view.at(0));
	if (mu != 0) {
	  int id = int(mu->globalTrack().id().id());
	  recLevelMap[id] = rec;
	}
      }
    }
  }
}

void RecLevelHelper::storeMatchMap(const edm::Event& event,
				   const string& mapName,
				   reco::CandViewMatchMap& map) const {
  edm::Handle<reco::CandViewMatchMap> matchMap;
  event.getByLabel(mapName, matchMap);

  if (matchMap.failedToGet()) {
    // Fail silently. The user will know what's wrong when absolutely
    // nothing gets matched.
    //edm::LogWarning("storeMatchMap")
    //  << "Couldn't get match map " << mapName << " from event!"
    //  << " Storing empty match map.";
    map = reco::CandViewMatchMap();
  }
  else
    map = *matchMap;
}

void RecLevelHelper::storeMatchMaps(const edm::Event& event) {
  edm::Handle<vector<int> > seeds;

  for (int i = 0; i < MAX_LEVELS; i++) {
    // Store MC truth match maps.
    storeMatchMap(event, "genMatch" + levelName(i), genMatchMap[i]);

    // Store the photon match maps (which only exist for global fit
    // leptons).
    if (i >= lGR)
      storeMatchMap(event, "photonMatch" + levelName(i), photonMatchMap[i]);
  }
}

int RecLevelHelper::id(const reco::CandidateRef& cand) const {
  if (cand.isNull())
    return -999;
  return cand.index();
}

int RecLevelHelper::id(const reco::CandidateBaseRef& cand) const {
  if (cand.isNull())
    return -999;
  return cand.key();
}

int RecLevelHelper::recLevel(const reco::CandidateBaseRef& cand) const {
  map<int,int>::const_iterator c = recLevelMap.find(cand.id().id());
  return c != recLevelMap.end() ? c->second : -1;
}

int RecLevelHelper::recLevel(const reco::CompositeCandidate& dil) const {
  int rec = -1;
  for (unsigned ilep = 0; ilep < dil.numberOfDaughters(); ilep++) {
    int r = recLevel(dileptonDaughter(dil, ilep));
    if (rec >= 0) {
      if (r != rec)
	return -1;
    }
    else
      rec = r;
  }
  return rec;
} 

int RecLevelHelper::originalRecLevel(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
  if (!isCocktailLevel(level))
    return level;

  // If we have a cocktail muon, then use its main track and the track
  // product id -> rec level map stored earlier to figure out what the
  // original rec level was.
  const reco::Muon* mu = toConcretePtr<reco::Muon>(cand);
  if (mu != 0) {
    int id = mu->globalTrack().id().id();
    map<int,int>::const_iterator c = recLevelMap.find(id);
    if (c != recLevelMap.end())
      return c->second;
  }

  // If we failed to find it, just return the cocktail level. (One way
  // to get get here is if cand is a cocktail muon produced by the
  // "official" method tevOptimized(), since it can return a GMR track
  // that has been refit once ("tevMuons_default"). Zprime2muAnalysis
  // does not currently use these tracks.)
  return level;
}

bool
RecLevelHelper::hasClosestPhoton(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
  checkRecLevel(level, "hasClosestPhoton");

  // No closest photon for non-offline fits.
  if (level < lGR)
    return false;

  const reco::CandViewMatchMap& mm = photonMatchMap[level];
  return mm.find(cand) != mm.end();
}

reco::Particle::LorentzVector
RecLevelHelper::closestPhoton(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
  checkRecLevel(level, "closestPhoton");

  // No closest photon for non-offline fits.
  if (level < lGR)
    return reco::Particle::LorentzVector();

  const reco::CandViewMatchMap& mm = photonMatchMap[level];
  // If no closest photon found, return a zero four-vector.
  if (mm.find(cand) == mm.end())
    return reco::Particle::LorentzVector();
  else
    return mm[cand]->p4();
}

const reco::CandidateBaseRef&
RecLevelHelper::matchGenLepton(const reco::CandidateBaseRef& lep) const {
  int oldlevel = recLevel(lep);
  checkRecLevel(oldlevel, "matchGenLepton");

  const int level = lGN;
  if (oldlevel == level)
    return lep;

  const reco::CandViewMatchMap& mm = genMatchMap[oldlevel];
  if (mm.find(lep) != mm.end())
    return mm[lep];
    
  return invalidRef;
}

const reco::CandidateBaseRef&
RecLevelHelper::matchOfflineLepton(const reco::CandidateBaseRef& lep,
				   const int level) const {
  checkRecLevel(level, "matchOfflineLepton");
  if (level < lGR)
    return invalidRef;

  /*
  // JMTBAD need better way to match different offline fits, this is
  // absolutely horrible.

  edm::View<reco::Candidate> leps;
  if (!get(*evt, level, lepInputs[level], leps))
    return invalidRef;

  const reco::Muon* orig = toConcretePtr<reco::Muon>(lep);
  if (orig == 0)
    return invalidRef;
  const reco::Track* orig_tk = orig->innerTrack().get();
  if (orig_tk == 0)
    return invalidRef;

  for (unsigned i = 0; i < leps.size(); ++i) {
    const reco::CandidateBaseRef& newlep = leps.refAt(i);

    const reco::Muon* mu = toConcretePtr<reco::Muon>(newlep);
    if (mu == 0)
      continue;
    const reco::Track* tk = mu->innerTrack().get();
    if (tk == 0)
      continue;

    if (tk->momentum() == orig_tk->momentum())
      return newlep;
  }
  */

  return invalidRef;
}
