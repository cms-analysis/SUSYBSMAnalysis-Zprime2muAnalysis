#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
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
  produces<CandidateCollection>();
}

void GenCandsFromSimTracksProducer::produce(Event& event,
					    const EventSetup& eSetup) {
  // Simulated tracks (i.e. GEANT particles).
  Handle<SimTrackContainer> simtracks;
  // and the associated vertices
  Handle<SimVertexContainer> simvertices;
  bool ok = true;
  try {
    event.getByLabel(src, simtracks);
    event.getByLabel(src, simvertices);
  } catch (...) {
    ok = false;
  }

  // make the output collection
  auto_ptr<CandidateCollection> cands(new CandidateCollection);

  if (ok && !simtracks.failedToGet() && !simvertices.failedToGet()) {
    // Add GEANT tracks which were not in PYTHIA list, since those in
    // the latter should already be in the GenParticleCandidates.
    for (SimTrackContainer::const_iterator isimtrk = simtracks->begin();
	 isimtrk != simtracks->end(); ++isimtrk) {
      int pythiaInd = isimtrk->genpartIndex();
      if (pythiaInd != -1) continue; // Skip PYTHIA tracks.

      const CLHEP::HepLorentzVector& pmu = isimtrk->momentum();

      // Make up a GenParticleCandidate from the GEANT track info.
      int charge = static_cast<int>(isimtrk->charge());
      Particle::LorentzVector p4;
      p4.SetXYZT(pmu.x(), pmu.y(), pmu.z(), pmu.t());
      Particle::Point vtx; // = (0,0,0) by default
      if (!isimtrk->noVertex()) {
	const CLHEP::HepLorentzVector& v
	  = (*simvertices)[isimtrk->vertIndex()].position();
	vtx = v;
      }
      int status = 1;
      GenParticleCandidate* genp
	= new GenParticleCandidate(charge, p4, vtx, isimtrk->type(),
				   status, true);
      cands->push_back(genp);
    }
  }
  else
    edm::LogWarning("GenCandsFromSimTracksProducer")
      << "no collection " << src << " in event; producing empty collection";

  event.put(cands);
}

DEFINE_FWK_MODULE(GenCandsFromSimTracksProducer);
