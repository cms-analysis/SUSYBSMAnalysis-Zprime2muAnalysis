#include "CommonTools/UtilAlgos/interface/Merger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

typedef Merger<reco::GenParticleCollection> GenParticleMerger;

DEFINE_FWK_MODULE(GenParticleMerger);
