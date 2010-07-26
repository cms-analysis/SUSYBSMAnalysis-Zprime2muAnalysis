import FWCore.ParameterSet.Config as cms

genSimLeptons = cms.EDProducer(
    "GenPlusSimParticleProducer",
    src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
    setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
    particleTypes = cms.vstring("mu+","e-"),
    filter        = cms.vstring('(abs(pdgId) == 13 || abs(pdgId) == 11) && pt > 2'),  # just for testing (optional)
    genParticles  = cms.InputTag("genParticles") # original genParticle list
    )

prunedGenSimLeptons = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticlesSimLep"),
    select = cms.vstring(
    "drop * ",
    "keep status <= 3",
    "++keep++ pdgId = {mu-} && pt > 2",
    "++keep++ pdgId = {e-} && pt > 2",
    #"+drop pdgId = {mu-}",
    #"+drop pdgId = {e-}"
    )
    )

allGenSimMergerSequence = cms.Sequence(genSimLeptons + prunedGenSimLeptons)

simParticles = cms.EDProducer(
    'PATGenCandsFromSimTracksProducer',
    src = cms.InputTag('g4SimHits'),
    setStatus = cms.int32(8),
    filter = cms.string('(abs(pdgId) == 13 || abs(pdgId) == 11) && pt > 2')
    )

genSimParticles = cms.EDProducer(
    'GenParticleMerger',
    src = cms.VInputTag(
    cms.InputTag('genParticles'),
    cms.InputTag('simParticles'))
    )

genSimMergerSequence = cms.Sequence(simParticles + genSimParticles)
