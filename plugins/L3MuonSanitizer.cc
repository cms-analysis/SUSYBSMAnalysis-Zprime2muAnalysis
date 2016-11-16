#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using namespace std;
using namespace edm;
using namespace reco;

class L3MuonSanitizer : public EDProducer {
public:
  explicit L3MuonSanitizer(const ParameterSet&);
  ~L3MuonSanitizer() {}

private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}

  InputTag src;
};

L3MuonSanitizer::L3MuonSanitizer(const ParameterSet& cfg) {
  src = cfg.getParameter<InputTag>("src");
  consumes<MuonTrackLinksCollection>(src);
  produces<MuonCollection>();
}

void L3MuonSanitizer::produce(Event& event,
				  const EventSetup& eSetup) {
  static const double muMass = 0.10566;

  // Try to get the L3 muons from the event.
  Handle<MuonTrackLinksCollection> muons;
  event.getByLabel(src, muons);

  // Make the output collection.
  auto_ptr<MuonCollection> cands(new MuonCollection);

  if (!muons.failedToGet()) {
    MuonTrackLinksCollection::const_iterator muon;
    for (muon = muons->begin(); muon != muons->end(); muon++) {
      // Null references to combined, tracker-only, and standalone muon tracks.
      TrackRef theTrack;
      TrackRef tkTrack;
      TrackRef muTrack;

      theTrack = muon->globalTrack();
      tkTrack  = muon->trackerTrack();
      muTrack  = muon->standAloneTrack();

      // Do the rest only if theTrack is not empty.
      if (theTrack.isNonnull()) {
	// make up a reco::Muon from the tracks
	Particle::Point vtx(theTrack->vx(), theTrack->vy(), theTrack->vz());
	Particle::LorentzVector p4;
	p4.SetXYZT(theTrack->px(), theTrack->py(), theTrack->pz(),
		   sqrt(theTrack->p()*theTrack->p() + muMass*muMass));
	Muon mu(theTrack->charge(), p4, vtx);
	mu.setGlobalTrack(theTrack);
	mu.setInnerTrack(tkTrack);
	mu.setOuterTrack(muTrack);
	cands->push_back(mu);
      }
    }
  }
  
  // Fail silently to prevent spamming messages. The user will know when
  // something is amiss from empty plots.

  //else
  //  edm::LogWarning("L3MuonSanitizer")
  //    << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

DEFINE_FWK_MODULE(L3MuonSanitizer);
