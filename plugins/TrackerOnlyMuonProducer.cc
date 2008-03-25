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

// JMTBAD this module should go away after understanding the
// difference between the output of option 1 of GlobalMuonProducer and
// just looking at the tracker tracks of these muons

class TrackerOnlyMuonProducer : public EDProducer {
public:
  explicit TrackerOnlyMuonProducer(const ParameterSet&);
  ~TrackerOnlyMuonProducer() {}

private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}

  InputTag src;
};

TrackerOnlyMuonProducer::TrackerOnlyMuonProducer(const ParameterSet& cfg) {
  src = cfg.getParameter<InputTag>("src");
  produces<MuonCollection>();
}

void TrackerOnlyMuonProducer::produce(Event& event, const EventSetup& eSetup) {
  // get the muons from the event
  Handle<MuonCollection> muons;
  bool ok = true;
  try {
    event.getByLabel(src, muons);
  } catch (...) {
    ok = false;
  }

  // make the output collection
  auto_ptr<MuonCollection> cands(new MuonCollection);
  const double MASS = 0.10566;

  if (ok && !muons.failedToGet()) {
    MuonCollection::const_iterator muon;
    for (muon = muons->begin(); muon != muons->end(); muon++) {
      // Null references to combined, tracker-only, and standalone muon tracks.
      TrackRef theTrack;
      TrackRef tkTrack;
      TrackRef muTrack;

      tkTrack = theTrack = muon->track();
      muTrack = muon->standAloneMuon();
      
      // Do the rest only if tkTrack is not empty.
      if (theTrack.isNonnull()) {
	// make up a real Muon from the tracks
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
    }
  }
  else
    edm::LogWarning("TrackerOnlyMuonProducer")
      << "no collection " << src << " in event; producing empty collection";
  
  event.put(cands);
}

DEFINE_FWK_MODULE(TrackerOnlyMuonProducer);
