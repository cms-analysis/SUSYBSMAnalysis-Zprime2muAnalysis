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

# Do some further pruning of the above-created genSimLeptons as
# detailed in the comments below.
prunedGenSimLeptons = cms.EDProducer('GenParticlePruner',
                                     src = cms.InputTag('genSimLeptons'),
                                     select = cms.vstring(
                                         'drop *',
                                         # Keep all doc-line particles (e.g. the Z0/Z', its parents, daughters, etc.).
                                         'keep status == 3',
                                         # Keep any final-state muon (either from the generator with status 1,
                                         # or from GEANT with status 8 as set above) that has potential to
                                         # reach the muon system (should reevaluate for muons produced far from
                                         # the IP), and all of their ancestors.
                                         '++keep abs(pdgId) == 13 && (status == 1 || status == 8)',
                                         # Keep any final-state non-GEANT electron and all of their ancestors.
                                         '++keep abs(pdgId) == 11 && status == 1',
                                         # Just to be sure, keep all the descendants of the Z0/Z'/G*.
                                         'keep++ pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039',
                                         )
                                     )

# A pruned list using only genParticles without merging the SIM-level
# particles. This is needed for AODSIM which doesn't have the latter.
prunedGenLeptons = prunedGenSimLeptons.clone(src = cms.InputTag('genParticles'))

genSequence = cms.Sequence(prunedGenLeptons)
genSimSequence = cms.Sequence(genSimLeptons + prunedGenSimLeptons)
