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
    edm::InputTag tag = config.getParameter<edm::InputTag>(tagname);
    if (tag.label() != "skip")
      lepInputs[i] = tag;
    else
      lepInputs[i] = edm::InputTag();

    tagname = "dileptons" + levelName(i, true);
    string diltag;
    diltag = config.getParameter<string>(tagname);
    if (diltag != "skip")
      dilInputs[i] = diltag;
    else
      dilInputs[i] = "";
  }
}

void RecLevelHelper::initEvent(const edm::Event& event) {
  for (int rec = 0; rec < MAX_LEVELS; rec++)
    warned[rec] = false;
  storeRecLevelMap(event);
  storeMatchMaps(event);
}

bool RecLevelHelper::getLeptonsView(const edm::Event& event, int level,
				    edm::View<reco::Candidate>& view) {
  edm::Handle<edm::View<reco::Candidate> > hview;
  event.getByLabel(lepInputs[level], hview);
  if (hview.failedToGet()) {
    // If there was no "best" collection in the event, try to use the
    // default global leptons (lgmr).
    if (level == lbest)
      return getLeptonsView(event, lgmr, view);
    
    if (!warned[level]) {
      string inp = lepInputs[level].encode();
      // Don't bother to warn about collections that are supposed to be missing.
      if (inp != ":" && inp != "") 
	edm::LogWarning("initEvent")
	  << "No event collection " << inp
	  << " found at rec level " << level << "; skipping";
      warned[level] = true;
    }
    return false;
  }

  view = *hview;
  return true;
}

bool RecLevelHelper::getLeptons(const edm::Event& event, int level,
				reco::CandidateBaseRefVector& leps,
				double ptMin) {
  edm::View<reco::Candidate> view;
  if (!getLeptonsView(event, level, view))
    return false;

  // Store refs to the leptons that are valid (i.e. pt > ptMin = 1e-3
  // by default).
  leps.clear();
  for (unsigned ilep = 0; ilep < view.size(); ilep++)
    if (view[ilep].pt() > ptMin)
      leps.push_back(view.refAt(ilep));

  return true;
}

bool RecLevelHelper::getDileptons(const edm::Event& event, int level, DilType type,
				  reco::CompositeCandidateCollection& coll) const {
  static const string extra[3] = {"", "Res", "Raw"};

  // No "resonances" for L1 or L2, so just get the plain dileptons.
  if (type == RES && (level == l1 || level == l2)) type = DIL;
  string label = dilInputs[level] + extra[type];

  edm::Handle<reco::CompositeCandidateCollection> hcoll;
  event.getByLabel(label, hcoll);
  if (hcoll.failedToGet()) {
    if (!warned[level]) {
      // Don't bother to warn about collections that are supposed to be missing.
      if (label != "")
	edm::LogWarning("initEvent")
	  << "No event collection " << label
	  << " found at rec level " << level << "; skipping";
    }
    return false;
  }

  coll = *hcoll;
  return true;
}

bool RecLevelHelper::recLevelOkay(const edm::Event& event, int level) {
  edm::View<reco::Candidate> view;
  return getLeptonsView(event, level, view);
}
  
string RecLevelHelper::makeMatchMapName(RecLevelHelper::MatchType mtype,
					const int irec,
					const int jrec) const {
  static const char* base[] = { "photonMatch", "closestMatch", "seedMatch" };
  string res = base[int(mtype)];
  res += levelName(irec, true);
  if (jrec >= 0)
    res += levelName(jrec, true);
  return res;
}

void RecLevelHelper::storeRecLevelMap(const edm::Event& event) {
  // We shouldn't have to re-store this every event...
  recLevelMap.clear();
  for (int rec = 0; rec < MAX_LEVELS; rec++) {
    edm::View<reco::Candidate> view;
    getLeptonsView(event, rec, view);

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
    edm::LogWarning("storeMatchMap")
      << "Couldn't get match map " << mapName << " from event!"
      << " Storing empty match map.";
    map = reco::CandViewMatchMap();
  }
  else
    map = *matchMap;
}

void RecLevelHelper::storeMatchMaps(const edm::Event& event) {
  edm::Handle<vector<int> > seeds;

  for (int i = 0; i < MAX_LEVELS; i++) {
    if (i >= l3)
      // Store the photon match maps (which only exist for global fit
      // leptons).
      storeMatchMap(event, makeMatchMapName(PHOTON, i), photonMatchMap[i]);

    // This inner loop does not include lbest, since we do not match
    // to that level, only from it.
    for (int j = 0; j < MAX_LEVELS-1; j++) {
      // Don't bother storing an identity map.
      if (i == j) continue;

      // Store closest match maps.
      storeMatchMap(event, makeMatchMapName(CLOSEST, i, j),
		    closestMatchMap[i][j]);

      // Also store by-seed match maps (which only exist for global
      // fits).
      if (i >= l3 && j >= l3)
	storeMatchMap(event, makeMatchMapName(SEED, i, j), seedMatchMap[i][j]);
    }
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
  for (level = ltk; level <= lpmr; level++) {
    // If the closest lepton at the other level has the same
    // four-vector, charge, and vertex, this level is from where the
    // cocktail muon came.
    const reco::CandidateBaseRef& clos = closestLepton(cand, level);
    if (clos.isNonnull() &&
	cand->p4() == clos->p4() && cand->charge() == clos->charge() &&
	cand->vertex() == clos->vertex())
      return level;
  }

  // If we didn't find an equal lepton in one of the above
  // collections, the lepton must have come from lgmr (either via
  // HEEPSelector or just copying the default GMR muons as the "best"
  // ones).
  return lgmr;
}

reco::Particle::LorentzVector
RecLevelHelper::closestPhoton(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
  checkRecLevel(level, "closestPhoton");
  // No closest photon for non-global fits.
  if (level < l3)
    return reco::Particle::LorentzVector();
  const reco::CandViewMatchMap& mm = photonMatchMap[level];
  // If no closest photon found, return a zero four-vector.
  if (mm.find(cand) == mm.end())
    return reco::Particle::LorentzVector();
  else
    return mm[cand]->p4();
}

const reco::CandidateBaseRef&
RecLevelHelper::matchLepton(const reco::CandidateBaseRef& lep,
			    const int level,
			    int whichMatch) const {
  int oldlevel = recLevel(lep);
  checkRecLevel(level, "matchLepton");
  checkRecLevel(oldlevel, "matchLepton");

  if (oldlevel == level)
    return lep;

  bool canMatchSeed = oldlevel >= l3 && level >= l3;
  bool doSeed;

  if (whichMatch == 0)
    doSeed = false;
  else if (whichMatch == 1) {
    if (!canMatchSeed) return invalidRef;
    doSeed = true;
  }
  else
    doSeed = canMatchSeed;

  const reco::CandViewMatchMap& mm = doSeed ? seedMatchMap[oldlevel][level]
                                         : closestMatchMap[oldlevel][level];
  if (mm.find(lep) != mm.end())
    return mm[lep];
  else if (whichMatch == -1 && doSeed) {
    const reco::CandViewMatchMap& cmm = closestMatchMap[oldlevel][level];
    if (cmm.find(lep) != cmm.end())
      return cmm[lep];
  }
    
  return invalidRef;
}
