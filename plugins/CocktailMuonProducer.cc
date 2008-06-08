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
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

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

  // Retrieve the matched-by-seed lepton from ltk to either lfms or
  // lpmr. (We only use a subset of the by-seed matching available, so
  // we just implement the tiny part we need here.)
  const CandidateBaseRef sameSeedLepton(const CandidateBaseRef& lep,
					const int level) const;

  // Get the chi2 prob for the muon's combined track.
  double cumulativeChiSquare(const CandidateBaseRef& mu) const;

  // Our implementation of cocktail methods for picking muons from the
  // various TeV muon reconstructors (either TMR, picking between
  // tracker-only and tracker+first muon station, or Piotr's, picking
  // between those two and also PMR); returns a reference to the one
  // picked.
  int cocktailMuon(const CandidateBaseRef& trk,
		   CandidateBaseRef& picked) const;

  // The driver routine which uses the cocktail method above to pick
  // "best" muons.
  bool findBestMuons(const Event& event,
		     const auto_ptr<MuonCollection>& bestMuons);

  bool debug;

  // Keep some statistics on what the cocktail picked.
  int ntrk, ngmr, nfms, npmr, ngpr, ntot;

  // Whether to use the TMR cocktail or Piotr's.
  bool useTMR;

  InputTag trackerOnlyMuons, toFMSMap, toPMRMap;
  Handle<CandViewMatchMap> toFMS, toPMR;
};

CocktailMuonProducer::CocktailMuonProducer(const ParameterSet& cfg)
  : debug(false),
    useTMR(cfg.getParameter<bool>("useTMR")),
    trackerOnlyMuons(cfg.getParameter<InputTag>("trackerOnlyMuons")),
    toFMSMap(cfg.getParameter<InputTag>("toFMSMap")),
    toPMRMap(cfg.getParameter<InputTag>("toPMRMap"))
{
  ntrk = ngmr = nfms = npmr = ngpr = ntot = 0;

  produces<MuonCollection>();
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

const CandidateBaseRef
CocktailMuonProducer::sameSeedLepton(const CandidateBaseRef& lep,
				     const int level) const {
  const Handle<CandViewMatchMap>& mm = level == lfms ? toFMS : toPMR;
  if (mm->find(lep) != mm->end())
    return (*mm)[lep];
  else
    return CandidateBaseRef(); // an invalid ref
}

double CocktailMuonProducer::cumulativeChiSquare(const CandidateBaseRef& mu) const {
  const TrackRef& track = toConcrete<Muon>(mu).combinedMuon();
  Genfun::CumulativeChiSquare myCumulativeChiSquare(int(track->ndof()));
  return -log(1.-myCumulativeChiSquare(track->chi2()));
}

int CocktailMuonProducer::cocktailMuon(const CandidateBaseRef& trk,
				       CandidateBaseRef& picked) const {
  const CandidateBaseRef& fms = sameSeedLepton(trk, lfms);
  const CandidateBaseRef& pmr = sameSeedLepton(trk, lpmr);
  int result = -1;

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

  if      (result == ltk)  picked = trk;
  else if (result == lfms) picked = fms;
  else if (result == lpmr) picked = pmr;
  
  return result;
}

bool CocktailMuonProducer::findBestMuons(const Event& event,
					 const auto_ptr<MuonCollection>& bestMuons) {
  // Start with tracker-only tracks.
  Handle<View<Candidate> > htkonly;
  try {
    event.getByLabel(trackerOnlyMuons, htkonly);
  } catch (const edm::Exception& e) {
  //  if (htkonly.failedToGet()) {
    edm::LogWarning("findBestMuons")
      << "Unable to get a view to the tracker-only muons; putting no"
      << " collection into event for 'best' muons.";
    return false;
  }

  if (debug) LogTrace("findBestMuons") << "Finding best muons for event #"
				       << event.id().event();

  for (unsigned imu = 0; imu < htkonly->size(); imu++) {
    const CandidateBaseRef& trkMu  = htkonly->refAt(imu);
    CandidateBaseRef bestMu;
    int picked = cocktailMuon(trkMu, bestMu);
    if (!bestMu.isNull())
      bestMuons->push_back(toConcrete<Muon>(bestMu));

    // Skip the statistics if we're not debugging.
    if (!debug) continue;

    if (picked >= 0) {
      if (picked == ltk) {
	ntrk++; LogTrace("findBestMuons") << " --> select Tracker only";
      }
      else if (picked == lgmr) {
	ngmr++; LogTrace("findBestMuons") << " --> select GMR";
	const CandidateBaseRef& pmrMu = sameSeedLepton(trkMu, lpmr);
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
  // Get the match maps. (Let the exception be thrown if the maps are
  // not found; this is a fatal error.)
  event.getByLabel(toFMSMap, toFMS);
  event.getByLabel(toPMRMap, toPMR);

  // Find the "best" muons according to the chosen cocktail, and put
  // them in the event.
  auto_ptr<MuonCollection> bestMuons(new MuonCollection);

  if (findBestMuons(event, bestMuons))
    event.put(bestMuons);
}

DEFINE_FWK_MODULE(CocktailMuonProducer);
