import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import HepPDTESSource

genSimLeptons = cms.EDProducer(
    "GenPlusSimParticleProducer",
    src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
    setStatus     = cms.int32(8),              # set status = 8 for GEANT GPs
    filter        = cms.vstring('(abs(pdgId) == 13 || abs(pdgId) == 11) && pt > 2'),
    genParticles  = cms.InputTag("genParticles") # original genParticle list
    )

prunedGenSimLeptons = cms.EDProducer('GenParticlePruner',
                                     src = cms.InputTag('genSimLeptons'),
                                     select = cms.vstring(
                                         'drop *',
                                         # keep all doc-line particles (e.g. the Z0/Z', its parents, daughters, etc.)
                                         'keep status == 3',
                                         # keep any final-state muon (either from the generator with status 1
                                         # or from GEANT with status 8 as set above) that has potential to
                                         # reach the muon system (should reevaluate for muons produced far from
                                         # the IP), and all of their ancestors
                                         '++keep abs(pdgId) == 13 && (status == 1 || status == 8)',
                                         # keep any final-state non-GEANT electron and all of their ancestors
                                         '++keep abs(pdgId) == 11 && status == 1',
                                         # just to be sure, keep all the descendants of the Z0/Z'/G*
                                         'keep++ pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039',
                                         )
                                     )

genSimSequence = cms.Sequence(genSimLeptons + prunedGenSimLeptons)
