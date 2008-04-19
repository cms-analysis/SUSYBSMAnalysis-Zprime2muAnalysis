#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Zprime2muAnalysis.h"

using namespace std;
using namespace edm;
using namespace reco;

class CocktailMuonProducer : public EDProducer {
public:
  explicit CocktailMuonProducer(const ParameterSet&);
  ~CocktailMuonProducer() {}

private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob();

  // Get the chi2 prob for the muon's combined track.
  double cumulativeChiSquare(const reco::CandidateBaseRef& mu) const;

  // Our implementation of cocktail methods for picking muons from the
  // various TeV muon reconstructors (either TMR, picking between
  // tracker-only and tracker+first muon station, or Piotr's, picking
  // between those two and also PMR); returns a reference to the one
  // picked.
  const reco::CandidateBaseRef&
    cocktailMuon(const reco::CandidateBaseRef& trk) const;

  // The driver routine which uses the cocktail method above to pick
  // "best" muons.
  bool findBestMuons(const Event& event,
		     const auto_ptr<MuonCollection>& bestMuons);

  bool debug;

  // Keep some statistics on what the cocktail picked.
  int ntrk, ngmr, nfms, npmr, ngpr, ntot;

  // Whether to use the TMR cocktail or Piotr's.
  bool useTMR;

  // If there are no extra TeV muon collections in the event, then
  // don't bother producing anything.
  bool useOtherMuonRecos;

  // Object to help with per-rec-level information.
  RecLevelHelper recLevelHelper;

  // An invalid CandidateBaseRef so that some of the methods above can
  // return by reference and still be able to return a null reference.
  reco::CandidateBaseRef invalidRef;

  // Keep track of which collection the cocktail chose. (Since we are
  // forced to copy the muons, downstream we will not be able to
  // determine which collections they originally came from, and
  // therefore will not be able to get their original rec level as
  // before.)
  vector<int> muonSources;
};

CocktailMuonProducer::CocktailMuonProducer(const ParameterSet& cfg)
  : debug(cfg.getUntrackedParameter<int>("verbosity", 0) > 0),
    useTMR(cfg.getParameter<bool>("useTMR")),
    useOtherMuonRecos(cfg.getParameter<bool>("useOtherMuonRecos"))
{
  recLevelHelper.init(cfg);
  ntrk = ngmr = nfms = npmr = ngpr = ntot = 0;

  produces<MuonCollection>();
  produces<vector<int> >();
}

void CocktailMuonProducer::endJob() {
  if (debug) {
    // Dump statistics on fractions of events from different
    // reconstructors picked by the cocktail.
    edm::LogInfo("CocktailMuonProducer")
      << "findBestMuons summary: "
      << "TO: "    << ntrk/double(ntot) << "; GMR: " << ngmr/double(ntot)
      << "; FMS: " << nfms/double(ntot) << "; PMR: " << npmr/double(ntot)
      << "; GMR=PMR: " << (ngmr > 0 ? ngpr/double(ngmr) : 0);
  }
}

double CocktailMuonProducer::cumulativeChiSquare(const reco::CandidateBaseRef& mu) const {
  const reco::TrackRef& track = toConcrete<reco::Muon>(mu).combinedMuon();
  Genfun::CumulativeChiSquare myCumulativeChiSquare(int(track->ndof()));
  return -log(1.-myCumulativeChiSquare(track->chi2()));
}


const reco::CandidateBaseRef&
CocktailMuonProducer::cocktailMuon(const reco::CandidateBaseRef& trk) const {
  const reco::CandidateBaseRef& fms = recLevelHelper.sameSeedLepton(trk, lfms);
  const reco::CandidateBaseRef& pmr = recLevelHelper.sameSeedLepton(trk, lpmr);
  int result = 0;

  bool trkValid = !trk.isNull();
  bool fmsValid = !fms.isNull();
  bool pmrValid = !pmr.isNull();

  double prob_trk = 1000.0, prob_fms = 1000.0, prob_pmr = 1000.0;

  if (trkValid) prob_trk = cumulativeChiSquare(trk);
  if (fmsValid) prob_fms = cumulativeChiSquare(fms);
  if (pmrValid) prob_pmr = cumulativeChiSquare(pmr);

  if (useTMR && trkValid && fmsValid) {
    // If we're doing the TMR cocktail, only choose between
    // tracker-only and tracker+first muon station.
    if (prob_fms - prob_trk > 30)
      result = ltk;
    else
      result = lfms;
  }
  else if (fmsValid && !pmrValid) {
    //if (trkValid && prob_fms - prob_trk > 30.) result = ltk;
    //else                                       result = lfms;
    result = lfms;
  }
  else if (!fmsValid && pmrValid) {
    result = lpmr;
    //if (trkValid && prob_pmr - prob_trk > 30.) result = ltk;
  }
  else if (fmsValid && pmrValid) {
    if (prob_fms - prob_pmr > 0.9) result = lpmr;
    //else {
    //  if (trkValid && prob_fms - prob_trk > 30.) result = ltk;
    //  else                                       result = lfms;
    //}
    else result = lfms;
  }

  if (debug)
    LogTrace("cocktailMuon")
      << "  Probabilities: trk = " << prob_trk
      << "; fms = " << prob_fms << "; pmr = " << prob_pmr
      << "; result = " << (useTMR ? "(TMR choice only) " : "") << result;

  if      (result == ltk)  return trk;
  else if (result == lfms) return fms;
  else if (result == lpmr) return pmr;
  else                     return invalidRef;
}

bool CocktailMuonProducer::findBestMuons(const Event& event,
					 const auto_ptr<MuonCollection>& bestMuons) {
  // Start with tracker-only tracks.
  View<Candidate> trackerOnlyMuons;
  if (!recLevelHelper.getView(event, ltk, trackerOnlyMuons)) {
    edm::LogWarning("findBestMuons")
      << "Unable to get a view to the tracker-only muons; putting no"
      << " collection into event for 'best' muons.";
    return false;
  }

  if (debug) LogTrace("findBestMuons") << "Finding best muons for event #"
				       << event.id().event();

  for (unsigned imu = 0; imu < trackerOnlyMuons.size(); imu++) {
    int picked = -1;
    const reco::CandidateBaseRef& trkMu = trackerOnlyMuons.refAt(imu);
    const reco::CandidateBaseRef& bestMu = cocktailMuon(trkMu);
    if (!bestMu.isNull()) {
      bestMuons->push_back(toConcrete<Muon>(bestMu));
      picked = recLevelHelper.recLevel(bestMu);
    }
    // Even if the muon was rejected and picked is still -1, keep that
    // information.
    muonSources.push_back(picked);

    if (!debug) continue;

    if (picked >= 0) {
      if (picked == ltk) {
	ntrk++; LogTrace("findBestMuons") << " --> select Tracker only";
      }
      else if (picked == lgmr) {
	ngmr++; LogTrace("findBestMuons") << " --> select GMR";
	const reco::CandidateBaseRef& pmrMu =
	  recLevelHelper.sameSeedLepton(trkMu, lpmr);
	if (!pmrMu.isNull() &&
	    fabs(bestMu->phi() - pmrMu->phi()) < 0.001 &&
	    fabs(bestMu->eta() - pmrMu->eta()) < 0.001 &&
	    fabs(bestMu->pt()  - pmrMu->pt()) < 0.001)
	  ngpr++;
      }
      else if (picked == lfms) {
	nfms++; LogTrace("findBestMuons") << " --> select TPFMS";
      }
      else if (picked == lpmr) {
	npmr++; LogTrace("findBestMuons") << " --> select PMR";
      }
      else {
	throw cms::Exception("findBestMuons") << "+++ Unknown outcome! +++\n";
      }
      ntot++;
    }
    else // i.e. an invalid result
      edm::LogWarning("findBestMuons") << " --> reject muon\n";
  }

  // Even if we rejected all muons, still put an empty collection into
  // the event so consumers know we were run.
  return true;
}

void CocktailMuonProducer::produce(Event& event,
				   const EventSetup& eSetup) {
  // If we're not to bother, put no collection into the event that way
  // consumers know to just get the default global muons.
  if (!useOtherMuonRecos)
    return;

  // Initialize things.
  recLevelHelper.initEvent(event);
  muonSources.clear();

  // Find the "best" muons according to the chosen cocktail, and put
  // them in the event.
  auto_ptr<MuonCollection> bestMuons(new MuonCollection);

  if (findBestMuons(event, bestMuons)) {
    event.put(bestMuons);

    // Store (by copying) original rec levels for later inspection.
    auto_ptr<vector<int> > sources(new vector<int>(muonSources));
    event.put(sources);
  }
}

DEFINE_FWK_MODULE(CocktailMuonProducer);
