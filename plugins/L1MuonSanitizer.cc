#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace l1extra;

// L1 takes phi from the edge of the bin, not the center.
// This needs to be corrected by 1.25 degrees or 0.0218 rad.
const double l1_phi_correction = 0.0218;
const double muMass = 0.10566;
const double ptMin = 1e-3;

class L1MuonSanitizer : public EDProducer {
public:
  explicit L1MuonSanitizer(const ParameterSet&);

private:
  virtual void produce(Event&, const EventSetup&);

  InputTag src;
};

L1MuonSanitizer::L1MuonSanitizer(const ParameterSet& cfg) {
  src = cfg.getParameter<InputTag>("src");
  produces<L1MuonParticleCollection>();
}

void L1MuonSanitizer::produce(Event& event,
				  const EventSetup& eSetup) {
  // Try to get the L1 muons from the event
  Handle<L1MuonParticleCollection> muons;
  event.getByLabel(src, muons);

  // make the output collection
  auto_ptr<L1MuonParticleCollection> cands(new L1MuonParticleCollection);

  if (!muons.failedToGet()) {
    l1extra::L1MuonParticleCollection::const_iterator muon;
    for (muon = muons->begin(); muon != muons->end(); muon++) {
      if (muon->pt() > ptMin) {
	math::PtEtaPhiMLorentzVector p4(muon->pt(),
					muon->eta(),
					muon->phi() + l1_phi_correction,
					muMass);
	math::XYZTLorentzVector cp4(p4.x(), p4.y(), p4.z(), p4.t());
	L1MuonParticle l1p(muon->charge(), cp4, muon->gmtMuonCand());
	cands->push_back(l1p);
      }
    }
  }
  else
    edm::LogWarning("L1MuonSanitizer")
      << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

DEFINE_FWK_MODULE(L1MuonSanitizer);
