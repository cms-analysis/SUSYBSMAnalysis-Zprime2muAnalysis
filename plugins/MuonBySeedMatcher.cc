#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using namespace std;
using namespace edm;
using namespace reco;

// Store same-seed matches for the input pair of globally-fit muon
// collections. The seeds are the stand-alone muon tracks. The matches
// are stored in the Event as a CandViewMatchMap. There's probably a
// better way to do this, since the collections for GMR, TK, FMS, PMR
// all seem sorted in the same order.

const int InvalidIndex = -999;

class MuonBySeedMatcher : public EDProducer {
public:
  explicit MuonBySeedMatcher(const ParameterSet&);

private:
  // Return the index into the seedTracks collection that this muon's
  // standAloneMuon (passed in as track) most closely matches.
  int matchMuonToSeed(const TrackRef& track,
		      const double maxMag2=1e99) const;

  // Perform by-seed matching from the muons in candsFrom to the ones
  // in candsTo and store the results in the pointed-to match map.
  void matchBySeed(const View<Candidate>& candsFrom,
		   const View<Candidate>& candsTo,
		   const auto_ptr<CandViewMatchMap>& matchMap) const;

  // The driver routine -- put the match maps into the event.
  virtual void produce(Event& event, const EventSetup& eSetup);

  bool debug;
  InputTag seedTracks, src, matched;
  Handle<TrackCollection> hseedTracks;
};

MuonBySeedMatcher::MuonBySeedMatcher(const ParameterSet& cfg)
  : debug(false),
    seedTracks(cfg.getParameter<InputTag>("seedTracks")),
    src(cfg.getParameter<InputTag>("src")),
    matched(cfg.getParameter<InputTag>("matched"))
{
  produces<CandViewMatchMap>();
}

int MuonBySeedMatcher::matchMuonToSeed(const TrackRef& track,
				       const double maxMag2) const {
  int closest = -1;
  double mag2, minMag2 = 1e99;
  
  ostringstream out;
  if (debug) out << "  Track to match:\n"
		 << "    charge = " << track->charge()
		 << " p = " << track->momentum() << endl
		 << "  Stand-alone muons:\n";
  
  TrackCollection::const_iterator seedmu = hseedTracks->begin();
  int ndx = 0;
  for ( ; seedmu != hseedTracks->end(); seedmu++, ndx++) {
    if (debug) out << "    charge = " << seedmu->charge()
		   << " p = " << seedmu->momentum();
    
    if (seedmu->charge() == track->charge()) {
      mag2 = (seedmu->momentum() - track->momentum()).mag2();
      if (debug) out << " mag2: " << mag2;
      
      if (mag2 < minMag2) {
	closest = ndx;
	minMag2 = mag2;
      }
    }
    
    if (debug) out << endl;
  }
  
  if (debug) {
    out << "  Closest is ndx = " << closest;
    LogVerbatim("matchMuonToSeed") << out.str();
  }

  if (minMag2 > maxMag2) {
    LogWarning("matchMuonToSeed")
      << "could not match stand-alone muon!\n"
      << "muon 3-momentum p = " << track->momentum();
    return InvalidIndex;
  }
  else
    return closest;
}

void MuonBySeedMatcher::matchBySeed(const View<Candidate>& candsFrom,
				    const View<Candidate>& candsTo,
				    const auto_ptr<CandViewMatchMap>& matchMap) const {
  // Match the two collections to their seed tracks, storing the index
  // into the track collection.
  vector<int> seedIndex[2];
  const View<Candidate>* cands[2] = { &candsFrom, &candsTo };
  for (unsigned i = 0; i < 2; i++) {
    for (unsigned ilep = 0; ilep < cands[i]->size(); ilep++) {
      const TrackRef& leptrk =
	(*cands[i])[ilep].get<TrackRef, StandAloneMuonTag>();
      if (leptrk.isNonnull())
	seedIndex[i].push_back(matchMuonToSeed(leptrk));
      else
	seedIndex[i].push_back(InvalidIndex);
    }
  }

  ostringstream out;
  if (debug) out << "Performing matchBySeed:\n";

  for (unsigned i = 0; i < candsFrom.size(); i++) {
    if (debug) out << "  from #" << setw(2) << i 
		   << " pt: " << candsFrom[i].pt()
		   << " eta: " << candsFrom[i].eta()
		   << " phi: " << candsFrom[i].phi()
		   << " seed: " << seedIndex[0][i] << endl;

    // Skip matching candidates with pt < PTMIN (defined in
    // Zprime2muAnalysis; perhaps we should cut those at config file
    // level?).
    if (candsFrom[i].pt() < 1e-3)
      continue;
    
    // Skip candidates that did not match to a seed.
    if (seedIndex[0][i] == InvalidIndex)
      continue;

    int foundIndex = -1;

    for (unsigned j = 0; j < candsTo.size(); j++) {
      if (debug) out << "    to #" << setw(2) << i 
		     << " pt: " << candsTo[j].pt()
		     << " eta: " << candsTo[j].eta()
		     << " phi: " << candsTo[j].phi()
		     << " seed: " << seedIndex[1][j] << endl;

      if (candsTo[j].pt() < 1e-3 ||
	  seedIndex[1][j] == InvalidIndex)
	continue;
      
      if (seedIndex[0][i] == seedIndex[1][j]) {
	// If two or more targets have the same seed, refuse to pick,
	// and let Zprime2muAnalysis default to picking by closest
	// match (since it will not find the entry in the seedMatchMap
	// it is looking for).
	if (foundIndex >= 0) {
	  foundIndex = -1;
	  break;
	}
	else
	  foundIndex = j;
      }
    }
    
    if (foundIndex >= 0) {
      if (debug) out << "   -> choosing " << foundIndex << endl;
      matchMap->insert(candsFrom.refAt(i), candsTo.refAt(foundIndex));
    }
  }

  if (debug) LogVerbatim("matchBySeed") << out.str();
}

void MuonBySeedMatcher::produce(Event& event,
				const EventSetup& eSetup) {
  // Get the seed tracks (the offline-fit stand-alone muons), and the
  // two muon collections.
  Handle<View<Candidate> > hsrc, hmatched;
  event.getByLabel(seedTracks, hseedTracks);
  event.getByLabel(src,        hsrc);
  event.getByLabel(matched,    hmatched);
  bool ok = !hseedTracks.failedToGet() && !hsrc.failedToGet()
    && !hmatched.failedToGet();

  if (debug) {
    ostringstream out;

    out << "Stand-alone tracks (offline muon system fit), used as seeds:\n";
    unsigned i = 0;
    for (TrackCollection::const_iterator seedmu = hseedTracks->begin();
	 seedmu != hseedTracks->end(); seedmu++, i++)
      out << "  #" << setw(3) << i << " charge: " << setw(2) << seedmu->charge()
	  << " p: " << seedmu->momentum() << endl;

    LogVerbatim("MuonBySeedMatcher") << out.str();
  }

  // Set up the output match map.
  auto_ptr<CandViewMatchMap> seedMatchMap(new CandViewMatchMap);

  if (ok)
    // Store by-seed matches for all globally-fit leptons.
    matchBySeed(*hsrc, *hmatched, seedMatchMap);
  else
    LogWarning("MuonBySeedMatcher")
      << "one of the input collections " << src << " or " << matched
      << " not in event; producing empty match map.";
  
  event.put(seedMatchMap);
}

DEFINE_FWK_MODULE(MuonBySeedMatcher);
