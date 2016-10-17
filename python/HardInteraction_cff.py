import FWCore.ParameterSet.Config as cms

hardInteraction = cms.PSet(src = cms.InputTag('prunedMCLeptons'),
                           doingElectrons = cms.bool(False),
                           allowFakeResonance = cms.bool(True),
                           resonanceIds = cms.vint32(32, 23, 39, 5000039),
                           shutUp = cms.bool(False),
                           )


hardInteraction_MiniAOD = cms.PSet(src = cms.InputTag('prunedGenParticles'),
                           doingElectrons = cms.bool(False),
                           allowFakeResonance = cms.bool(True),
                           resonanceIds = cms.vint32(32, 23, 39, 5000039),
                           shutUp = cms.bool(False),
                           )
