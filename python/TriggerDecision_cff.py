import FWCore.ParameterSet.Config as cms

triggerDecision = cms.PSet(doingElectrons = cms.bool(False),
                           useTrigger = cms.bool(True),
                           l1GtObjectMap = cms.InputTag('hltL1GtObjectMap'),
                           hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
                           l1Paths = cms.vstring('L1_SingleMu12'),
                           hltPaths = cms.vstring('HLT_Mu30_v1')
                           )
