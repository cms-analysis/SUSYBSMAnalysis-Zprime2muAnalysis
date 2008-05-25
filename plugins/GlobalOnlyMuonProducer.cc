#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using namespace std;
using namespace edm;
using namespace reco;

// Module to copy only the GlobalMuons out of the default
// MuonCollection (and ignore the Tracker, CaloMuons).

class GlobalOnlyMuonProducer : public EDProducer {
public:
  explicit GlobalOnlyMuonProducer(const ParameterSet&);
  ~GlobalOnlyMuonProducer() {}

private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}

  InputTag src;

  // Allow building the muon from just the tracker track. This
  // functionality should go away after understanding the difference
  // between the output of option 1 of GlobalMuonProducer and just
  // looking at the tracker tracks of these muons
  bool fromTrackerTrack;
};

GlobalOnlyMuonProducer::GlobalOnlyMuonProducer(const ParameterSet& cfg)
  : src(cfg.getParameter<InputTag>("src")),
    fromTrackerTrack(cfg.getParameter<bool>("fromTrackerTrack"))
{
  produces<MuonCollection>();
}

void GlobalOnlyMuonProducer::produce(Event& event, const EventSetup& eSetup) {
  // get the muons from the event
  Handle<MuonCollection> muons;
  event.getByLabel(src, muons);

  // make the output collection
  auto_ptr<MuonCollection> cands(new MuonCollection);
  static const double MASS = 0.10566;

  if (!muons.failedToGet()) {
    MuonCollection::const_iterator muon;
    for (muon = muons->begin(); muon != muons->end(); muon++) {
      if (muon->isGlobalMuon()) {
	if (fromTrackerTrack) {
	  TrackRef theTrack = muon->track();
	  TrackRef tkTrack = muon->track();
	  TrackRef muTrack = muon->standAloneMuon();
	  
	  // Make up a real Muon from the tracker track.
	  Particle::Point vtx(theTrack->vx(), theTrack->vy(), theTrack->vz());
	  Particle::LorentzVector p4;
	  p4.SetXYZT(theTrack->px(), theTrack->py(), theTrack->pz(),
		     sqrt(theTrack->p()*theTrack->p() + MASS*MASS));
	  //Muon mu(theTrack->charge(), p4, vtx);
	  Muon* mu = muon->clone();
	  mu->setCharge(theTrack->charge());
	  mu->setP4(p4);
	  mu->setVertex(vtx);
	  mu->setCombined(theTrack);
	  mu->setTrack(tkTrack);
	  mu->setStandAlone(muTrack);
	  cands->push_back(*mu);
	}
	else
	  cands->push_back(*muon);
      }
    }
  }
  else
    edm::LogWarning("GlobalOnlyMuonProducer")
      << "no collection " << src << " in event; producing empty collection";
  
  event.put(cands);
}

DEFINE_FWK_MODULE(GlobalOnlyMuonProducer);
