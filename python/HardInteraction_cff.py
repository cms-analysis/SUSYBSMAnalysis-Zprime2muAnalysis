import FWCore.ParameterSet.Config as cms

hardInteraction = cms.PSet(src = cms.InputTag('genParticles'),
                           doingElectrons = cms.bool(False),
                           allowFakeResonance = cms.bool(True),
                           resonanceIds = cms.vint32(32, 23, 39, 5000039),
                           )
