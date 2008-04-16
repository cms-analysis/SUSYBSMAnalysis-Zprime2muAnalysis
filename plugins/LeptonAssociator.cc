#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

using namespace std;
using namespace edm;
using namespace reco;

// Store closest (in delta R) matches for all muons at all rec levels,
// store same-seed matches for all muons and all pairs of different
// rec levels above and including L3 (i.e. all global fits), and store
// matches for all global fits to reconstructed photons. These matches
// are stored in the Event as CandViewMatchMaps.

class LeptonAssociator : public EDProducer {
public:
  explicit LeptonAssociator(const ParameterSet&);

private:
  // Match each of the candsFrom to the closest one of candsTo, where
  // the metric is in eta-phi space (delta R). We do something similar
  // to TrivialDeltaRMatcher, only do it from code so that it can be
  // called in our loop. Require either that delta R is less than the
  // specified value (i.e. within a circle in eta-phi space), or that
  // each of delta eta and delta phi are less than the specified value
  // (i.e. within a rectangle). The latter option is a throwback to
  // old code, and is only used when matching leptons to leptons.
  // jrec == -1 signifies that the "to" collection is photons.
  void minDeltaRMatch(const auto_ptr<CandViewMatchMap>& matchMap,
		      const View<Candidate>& candsFrom,
		      const View<Candidate>& candsTo,
		      const int irec, const int jrec) const;

  // Return the index into the standAloneMuons collection staTracks
  // that this Muon's standAloneMuon (passed in as track) most closely
  // matches -- use to mock up seedIndex.
  int matchStandAloneMuon(const TrackRef& track,
			  bool relaxedMatch=true) const;

  // For all global fit muons (i.e. L3 and above), find the "seed"
  // track in the standAloneMuons collection, and store the index to
  // it in the seedIndex vector. This will be used in matchBySeed()
  // below.
  void matchLeptonsToSeeds();

  // Perform by-seed matching from leptons at rec level irec to those
  // at level jrec, using the seed indices found already, and store
  // the results in the pointed-to match map.
  void matchBySeed(const auto_ptr<CandViewMatchMap>& matchMap,
		   const int irec, const int jrec) const;

  // Clear out any containers and get the leptons, photons, and seed
  // tracks from the event.
  void initialize(const Event& event);

  // The driver routine -- put the match maps into the event.
  virtual void produce(Event&, const EventSetup&);

  bool debug;
  bool doingElectrons;
  InputTag standAloneMuons;
  InputTag photonInput;

  View<Candidate> leptons[MAX_LEVELS];
  View<Candidate> photons;
  TrackCollection seedTracks;
  vector<int> seedIndex[MAX_LEVELS];

  // Object to help with per-rec-level information.
  RecLevelHelper recLevelHelper;
};

LeptonAssociator::LeptonAssociator(const ParameterSet& cfg)
  : debug(false),
    doingElectrons(cfg.getParameter<bool>("doingElectrons")),
    standAloneMuons(cfg.getParameter<InputTag>("standAloneMuons")),
    photonInput(cfg.getParameter<InputTag>("photons"))
{
  recLevelHelper.init(cfg);

  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    if (irec >= l3) {
      produces<CandViewMatchMap>(recLevelHelper.makeMatchMapName(RecLevelHelper::PHOTON, irec));
      produces<vector<int> >(recLevelHelper.makeMatchMapName(RecLevelHelper::SEED, irec));
    }
    for (int jrec = 0; jrec < MAX_LEVELS; jrec++) {
      produces<CandViewMatchMap>(recLevelHelper.makeMatchMapName(RecLevelHelper::CLOSEST, irec, jrec));
      if (irec >= l3 && jrec >= l3)
	produces<CandViewMatchMap>(recLevelHelper.makeMatchMapName(RecLevelHelper::SEED, irec, jrec));
    }
  }
}

void
LeptonAssociator::minDeltaRMatch(const auto_ptr<CandViewMatchMap>& matchMap,
				 const View<Candidate>& candsFrom,
				 const View<Candidate>& candsTo,
				 const int irec, const int jrec) const {
  const bool leptonCuts = jrec >= 0;
  const double maxMinDist = leptonCuts ? 0.5 : 999;
  double dphi, deta, dist;

  ostringstream out;
  if (debug) {
    out << "Performing minDeltaRMatch " << irec << " -> ";
    if (jrec < 0) out << "photons";
    else out << jrec;
    out << ":\n";
  }

  for (unsigned i = 0; i < candsFrom.size(); i++) {
    if (debug) out << "  from #" << setw(2) << i 
		   << " pt: " << candsFrom[i].pt()
		   << " eta: " << candsFrom[i].eta()
		   << " phi: " << candsFrom[i].phi() << endl;
    // Skip matching candidates with pt < PTMIN (defined in
    // Zprime2muAnalysis; perhaps we should cut those at config file
    // level?).
    if (candsFrom[i].pt() < 1e-3)
      continue;

    double minDist = 999;
    int closest = -1;
    for (unsigned j = 0; j < candsTo.size(); j++) {
      if (debug) out << "    to #" << setw(2) << j
		     << " pt: " << candsTo[j].pt()
		     << " eta: " << candsTo[j].eta()
		     << " phi: " << candsTo[j].phi();

      double jpt = candsTo[j].pt();
      // Also do not match to generator-level leptons with too low pt.
      if (leptonCuts && ((jrec == lgen && jpt < 1) || jpt < 1e-3)) {
	if (debug) out << endl;
	continue;
      }

      // Calculate dphi and deta ourselves so the rectangle cut is
      // possible.
      dphi = fabs(deltaPhi(candsFrom[i].phi(), candsTo[j].phi()));
      deta = fabs(candsFrom[i].eta() - candsTo[j].eta());
      // dist = deltaR(candsFrom[i], candsTo[j]);
      dist = sqrt(dphi*dphi + deta*deta);

      bool passCut = leptonCuts ? dphi < maxMinDist && deta < maxMinDist
	: dist < maxMinDist;
      if (passCut && dist < minDist) {
	minDist = dist;
	closest = int(j);
      }
      
      if (debug) out << " dist: " << dist << " minDist: " << minDist
		     << " closest: " << closest << endl;
    }
    
    if (closest >= 0) {
      if (debug) out << "   -> choosing " << closest << endl;
      matchMap->insert(candsFrom.refAt(i), candsTo.refAt(closest));
    }
    else {
      if (debug) out << "   -> did not find a match!\n";
    }
  }
  
  if (debug) LogVerbatim("minDeltaRMatch") << out.str();
}

int LeptonAssociator::matchStandAloneMuon(const TrackRef& track,
					  bool relaxedMatch) const {
  // JMTBAD there's probably a better way to do this, since the
  // collections for GMR, TK, FMS, PMR all seem sorted in the same way
  int closest = -1;
  double mag2, minMag2 = 1e99;
  const double maxMag2 = 1; // more than enough for a comparison
  // between doubles
  int ndx = 0;
  ostringstream out;
  if (debug) out << "  Track to match:\n"
		 << "    charge = " << track->charge()
		 << " p = " << track->momentum() << endl
		 << "  Stand-alone muons:\n";
  for (TrackCollection::const_iterator seedmu = seedTracks.begin();
       seedmu != seedTracks.end(); seedmu++) {
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
    ndx++;
  }

  if (debug) out << "  Closest is ndx = " << closest;

  if (debug) LogVerbatim("matchStandAloneMuon") << out.str();

  if (minMag2 > maxMag2 && !relaxedMatch) {
    LogWarning("matchStandAloneMuon")
      << "could not match stand-alone muon!\n"
      << "muon 3-momentum p = " << track->momentum();
    return -999;
  }
  else
    return closest;
}

void LeptonAssociator::matchLeptonsToSeeds() {
  for (int irec = l3; irec < MAX_LEVELS; irec++) {
    for (unsigned ilep = 0; ilep < leptons[irec].size(); ilep++) {
      const RecoCandidate& lep
      	= toConcrete<RecoCandidate>(leptons[irec][ilep]);
      const TrackRef& leptrk = lep.standAloneMuon();
      if (debug) LogVerbatim("matchLeptonsToSeeds")
	<< "Matching to stand-alone muons, rec level " << irec
	<< " lepton #" << ilep << ":";
      int ndx = matchStandAloneMuon(leptrk);
      seedIndex[irec].push_back(ndx);
    }
  }
}
    
void LeptonAssociator::matchBySeed(const auto_ptr<CandViewMatchMap>& matchMap,
				   const int irec, const int jrec) const {
  if (irec < l3 || jrec < l3)
    throw cms::Exception("matchBySeed") << "Cannot match " << irec
					<< " -> " << jrec << " by seed!\n";
      
  const View<Candidate>& candsFrom = leptons[irec];
  const View<Candidate>& candsTo = leptons[jrec];

  ostringstream out;
  if (debug) out << "Performing matchBySeed " << irec << " -> "
		 << jrec << ":\n";

  for (unsigned i = 0; i < candsFrom.size(); i++) {
    if (debug) out << "  from #" << setw(2) << i 
		   << " pt: " << candsFrom[i].pt()
		   << " eta: " << candsFrom[i].eta()
		   << " phi: " << candsFrom[i].phi()
		   << " seed: " << seedIndex[irec][i] << endl;

    // Skip matching candidates with pt < PTMIN (defined in
    // Zprime2muAnalysis; perhaps we should cut those at config file
    // level?).
    if (candsFrom[i].pt() < 1e-3)
      continue;
    
    // Skip candidates that did not match to a seed.
    if (seedIndex[irec][i] == -999)
      continue;

    int foundIndex = -1;

    for (unsigned j = 0; j < candsTo.size(); j++) {
      if (debug) out << "    to #" << setw(2) << i 
		     << " pt: " << candsTo[j].pt()
		     << " eta: " << candsTo[j].eta()
		     << " phi: " << candsTo[j].phi()
		     << " seed: " << seedIndex[jrec][j] << endl;

      if (candsTo[j].pt() < 1e-3 ||
	  seedIndex[jrec][j] == -999)
	continue;
      
      if (seedIndex[irec][i] == seedIndex[jrec][j]) {
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

void LeptonAssociator::initialize(const Event& event) {
  // We're storing views by value, so keep an empty one handy to copy
  // from in case the below getByLabels fail (else we end up somehow
  // using the view from the last event!).
  View<Candidate> emptyView;

  // Clear out seedIndex.
  for (int irec = 0; irec < MAX_LEVELS; irec++)
    seedIndex[irec].clear();

  // Store views to the leptons from the event.
  Handle<View<Candidate> > hview;
  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    View<Candidate> view;
    bool ok = recLevelHelper.getView(event, irec, view);
    leptons[irec] = ok ? view : emptyView;
  }

  // Store a view to the photons as well.
  bool ok = true;
  try {
    event.getByLabel(photonInput, hview);
  } catch (const cms::Exception& e) {
    //if (hview.failedToGet()) {
    ok = false;
  }
  photons = ok ? *hview : emptyView;

  // Store the seed tracks (the offline-fit stand-alone muons).
  Handle<TrackCollection> hcoll;
  event.getByLabel(standAloneMuons, hcoll);
  seedTracks = *hcoll;
}
				   
void LeptonAssociator::produce(Event& event,
			       const EventSetup& eSetup) {
  if (debug)
    LogInfo("LeptonAssociator") << "Producing match maps and 'best' leptons!";

  // Set up objects we're going to look at.
  initialize(event);
  
  // Match all global fits to their seed tracks. (Not currently
  // implemented for electrons since there is only one global
  // electron.)
  if (!doingElectrons) matchLeptonsToSeeds();

  for (int irec = 0; irec < MAX_LEVELS; irec++) {
    if (irec >= l3) {
      // Store closest photon for all global fit leptons.
      auto_ptr<CandViewMatchMap> photonMatchMap(new CandViewMatchMap);
      minDeltaRMatch(photonMatchMap, leptons[irec], photons,
		     irec, -1);
      event.put(photonMatchMap,
		recLevelHelper.makeMatchMapName(RecLevelHelper::PHOTON, irec));

      // Store (by copying) seed indices for later inspection.
      auto_ptr<vector<int> > seeds(new vector<int>(seedIndex[irec]));
      event.put(seeds,
		recLevelHelper.makeMatchMapName(RecLevelHelper::SEED, irec));
    }

    for (int jrec = 0; jrec < MAX_LEVELS; jrec++) {
      if (irec != jrec) {
	// Store closest matches for all pairs of rec levels.
	auto_ptr<CandViewMatchMap> closestMatchMap(new CandViewMatchMap);
	minDeltaRMatch(closestMatchMap, leptons[irec], leptons[jrec],
		       irec, jrec);
	event.put(closestMatchMap,
		  recLevelHelper.makeMatchMapName(RecLevelHelper::CLOSEST,
						  irec, jrec));

	if (irec >= l3 && jrec >= l3) {
	  // Store by-seed matches for all globally-fit leptons.
	  auto_ptr<CandViewMatchMap> seedMatchMap(new CandViewMatchMap);
	  matchBySeed(seedMatchMap, irec, jrec);
	  event.put(seedMatchMap,
		    recLevelHelper.makeMatchMapName(RecLevelHelper::SEED,
						    irec, jrec));
	}
      }
    }
  }
}

DEFINE_FWK_MODULE(LeptonAssociator);
