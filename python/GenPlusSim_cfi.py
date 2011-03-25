import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import HepPDTESSource

# Take the SIM-level GEANT-produced particles in the g4SimHits
# collections (or famosSimHits from FastSim) and create
# reco::GenParticles out of them for transparent use (their statuses
# set to 8 so they can be distinguished easily), and merge these with
# the PYTHIA-produced (or whatever GEN level particles) genParticles
# to make a new collection. Only keep electrons and muons.
genSimLeptons = cms.EDProducer(
    'GenPlusSimParticleProducer',
    src           = cms.InputTag('g4SimHits'), # use 'famosSimHits' for FastSim
    setStatus     = cms.int32(8),              # set status = 8 for GEANT GPs
    filter        = cms.vstring('(abs(pdgId) == 13 || abs(pdgId) == 11) && pt > 2'), # keep only muons and electrons with pT greater than 2 GeV
    genParticles  = cms.InputTag('genParticles') # original genParticle list
    )
