#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"

using namespace std;

void RecLevelHelper::init(const edm::ParameterSet& config) {
  doingElectrons = config.getParameter<bool>("doingElectrons");

  if (doingElectrons) {
    // JMTBAD preserve the order until L1/HLT electrons are implemented
    inputs = config.getParameter<vector<edm::InputTag> >("elInputs");
    while (int(inputs.size()) < MAX_LEVELS)
      inputs.push_back(edm::InputTag());
    inputs[lgmr] = inputs[l1];
    inputs[l1] = edm::InputTag();
  }
  else
    inputs = config.getParameter<vector<edm::InputTag> >("muInputs");

  if (int(inputs.size()) < MAX_LEVELS) {
    edm::LogWarning("") << "inputs.size() < MAX_LEVELS; padding it";
    while (int(inputs.size()) < MAX_LEVELS)
      inputs.push_back(edm::InputTag());
  }
}

void RecLevelHelper::initEvent(const edm::Event& event) {
  for (int rec = 0; rec < MAX_LEVELS; rec++)
    warned[rec] = false;
  storeRecLevelMap(event);
  storeMatchMaps(event);
}

bool RecLevelHelper::getView(const edm::Event& event,
			     int level,
			     edm::View<reco::Candidate>& view) {
  edm::Handle<edm::View<reco::Candidate> > hview;
  try {
    event.getByLabel(inputs[level], hview);
  } catch (const cms::Exception& e) {
  //if (hview.failedToGet()) {
    if (!warned[level]) {
      edm::LogWarning("initEvent")
	<< "No event collection " << inputs[level]
	<< " found at rec level " << level << "; skipping";
      warned[level] = true;
    }
    return false;
  }
  view = *hview;
  return true;
}

string RecLevelHelper::makeMatchMapName(RecLevelHelper::MatchType mtype,
					const int irec,
					const int jrec) const {
  static const char* base[] = { "photon", "closest", "seed" };
  string res = base[int(mtype)];
  res += str_level_short[irec];
  if (jrec >= 0)
    res += str_level_short[jrec];
  return res;
}

void RecLevelHelper::storeRecLevelMap(const edm::Event& event) {
  // We shouldn't have to re-store this every event...
  recLevelMap.clear();
  for (int rec = 0; rec < MAX_LEVELS; rec++) {
    edm::View<reco::Candidate> view;
    getView(event, rec, view);

    // Cache the product ids for each collection so we can look up the
    // rec level from the CandidateBaseRef.
    int key = int(view.id().id());
    if (recLevelMap.find(key) == recLevelMap.end())
      recLevelMap[key] = rec;
  }
}

void RecLevelHelper::storeMatchMaps(const edm::Event& event) {
  edm::Handle<reco::CandViewMatchMap> matchMap;
  edm::Handle<vector<int> > seeds;
  for (int i = 0; i < MAX_LEVELS; i++) {
    if (i >= l3) {
      // Store the photon match maps (which only exist for global fit
      // leptons).
      event.getByLabel("leptonMatches", makeMatchMapName(PHOTON, i), matchMap);
      photonMatchMap[i] = *matchMap;

      // Store the seed indices, which map candidates to their seeds
      // via the index into the seedIndices vector.
      event.getByLabel("leptonMatches", makeMatchMapName(SEED, i), seeds);
      seedIndices[i] = *seeds;
    }

    for (int j = 0; j < MAX_LEVELS; j++) {
      // Don't bother storing an identity map.
      if (i == j) continue;

      // Store closest and by-seed match maps (the latter only
      // existing for global fits).
      string sublabel = str_level_short[i] + str_level_short[j];
      event.getByLabel("leptonMatches", makeMatchMapName(CLOSEST, i, j),
		       matchMap);
      closestMatchMap[i][j] = *matchMap;
      if (i >= l3 && j >= l3) {
	event.getByLabel("leptonMatches", makeMatchMapName(SEED, i, j),
			 matchMap);
	seedMatchMap[i][j] = *matchMap;
      }
    }
  }
}

int RecLevelHelper::recLevel(const reco::CandidateBaseRef& cand) const {
  map<int,int>::const_iterator c = recLevelMap.find(cand.id().id());
  if (c != recLevelMap.end()) return c->second;
  else return -1;
}

void
RecLevelHelper::checkRecLevel(const int level, const char* name) const {
  if (level < 0 || level > MAX_LEVELS)
    throw cms::Exception(name)
      << "invalid level " << level << " in " << name << "!\n";
}

reco::Particle::LorentzVector
RecLevelHelper::closestPhoton(const reco::CandidateBaseRef& cand) const {
  int level = recLevel(cand);
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

int RecLevelHelper::seedIndex(const reco::CandidateBaseRef& cand) const {
  return seedIndices[recLevel(cand)][cand.key()];
}

const reco::CandidateBaseRef&
RecLevelHelper::matchLepton(const reco::CandidateBaseRef& lep,
			       const int level,
			       int whichMatch) const {
  int oldlevel = recLevel(lep);

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

const reco::CandidateBaseRef&
RecLevelHelper::closestLepton(const reco::CandidateBaseRef& lep,
				 const int level) const {
  return matchLepton(lep, level, 0);
}

const reco::CandidateBaseRef&
RecLevelHelper::sameSeedLepton(const reco::CandidateBaseRef& lep,
				  const int level) const {
  return matchLepton(lep, level, 1);
}

const reco::CandidateBaseRef&
RecLevelHelper::matchedLepton(const reco::CandidateBaseRef& lep,
				 const int level) const {
  return matchLepton(lep, level, -1);
}
