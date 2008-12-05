#include <string>

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"

using namespace std;

void checkRecLevel(const int level, const char* name) {
  if (level < 0 || level > MAX_LEVELS)
    throw cms::Exception(name)
      << "invalid level " << level << " in " << name << "!\n";
}

const std::string& levelName(const int rec, bool shortVersion) {
  checkRecLevel(rec, "levelName");
  if (shortVersion) return levelNamesShort[rec];
  else return levelNames[rec];
}

void RecLevelHelper::init(const edm::ParameterSet& config) {
  for (int i = 0; i < MAX_LEVELS; i++) {
    string tagname = "leptons" + levelName(i, true);
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
    // If there was no "best" collection in the event, try to use the
    // default global leptons (lgmr).
    if (level == lbest)
      return get(event, lgmr, tag, coll);
    
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
    storeMatchMap(event, "genMatch" + levelName(i, true), genMatchMap[i]);

    // Store the photon match maps (which only exist for global fit
    // leptons).
    if (i >= lgmr)
      storeMatchMap(event, "photonMatch" + levelName(i, true), photonMatchMap[i]);
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
	return lbest;
    }
    else
      rec = r;
  }
  return rec;
}  

int RecLevelHelper::originalRecLevel(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
  if (level != lbest)
    return level;

  // Try to look in the original collections the cocktail muon chose
  // from to see which level it came from.
  for (level = lgmr; level <= lpmr; level++) {
    // If the closest lepton at the other level has the same
    // four-vector, charge, and vertex, this level is from where the
    // cocktail muon came.
    const reco::CandidateBaseRef& clos = matchOfflineLepton(cand, level);
    if (clos.isNonnull() &&
	cand->p4() == clos->p4() && cand->charge() == clos->charge() &&
	cand->vertex() == clos->vertex())
      return level;
  }
  
  return lbest;
}

reco::Particle::LorentzVector
RecLevelHelper::closestPhoton(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
  checkRecLevel(level, "closestPhoton");

  // No closest photon for non-offline fits.
  if (level < lgmr)
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

  const int level = lgen;
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
  if (level < lgmr)
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
