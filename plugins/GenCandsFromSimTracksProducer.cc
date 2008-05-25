#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

using namespace std;
using namespace edm;
using namespace reco;

class GenCandsFromSimTracksProducer : public EDProducer {
public:
  explicit GenCandsFromSimTracksProducer(const ParameterSet&);
  ~GenCandsFromSimTracksProducer() {}

private:
  virtual void beginJob(const EventSetup&) {}
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() {}

  InputTag src;
};

GenCandsFromSimTracksProducer::GenCandsFromSimTracksProducer(
					        const ParameterSet& cfg) {
  src = cfg.getParameter<InputTag>("src");
  produces<GenParticleCollection>();
}

void GenCandsFromSimTracksProducer::produce(Event& event,
					    const EventSetup& eSetup) {
  // Simulated tracks (i.e. GEANT particles).
  Handle<SimTrackContainer> simtracks;
  event.getByLabel(src, simtracks);

  // and the associated vertices
  Handle<SimVertexContainer> simvertices;
  event.getByLabel(src, simvertices);

  // make the output collection
  auto_ptr<GenParticleCollection> cands(new GenParticleCollection);

  if (!simtracks.failedToGet() && !simvertices.failedToGet()) {
    // Add GEANT tracks which were not in PYTHIA list, since those in
    // the latter should already be in the GenParticleCandidates.
    for (SimTrackContainer::const_iterator isimtrk = simtracks->begin();
	 isimtrk != simtracks->end(); ++isimtrk) {
      int pythiaInd = isimtrk->genpartIndex();
      if (pythiaInd != -1) continue; // Skip PYTHIA tracks.

      // Make up a GenParticleCandidate from the GEANT track info.
      int charge = static_cast<int>(isimtrk->charge());
      Particle::LorentzVector p4 = isimtrk->momentum();
      Particle::Point vtx; // = (0,0,0) by default
      if (!isimtrk->noVertex())
	vtx = (*simvertices)[isimtrk->vertIndex()].position();
      int status = 1;

      GenParticle genp(charge, p4, vtx, isimtrk->type(), status, true);
      cands->push_back(genp);
    }
  }
  else
    edm::LogWarning("GenCandsFromSimTracksProducer")
      << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

DEFINE_FWK_MODULE(GenCandsFromSimTracksProducer);
