import FWCore.ParameterSet.Config as cms

muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRDPtLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    andOr                 = cms.bool(False),
    filterIdsEnum         = cms.vstring('TriggerMuon'),
    filterIds             = cms.vint32(0),
    filterLabels          = cms.vstring('*'),
    pathNames             = cms.vstring('HLT_Mu9', 'HLT_Mu11', 'HLT_DoubleMu3', 'HLT_DoubleMu3_v2', 'HLT_Mu13_v1', 'HLT_Mu15_v1'),
    collectionTags        = cms.vstring('*'),
    # Follow VBTF's matching criteria.
    maxDPtRel             = cms.double(1),
    maxDeltaR             = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)
)
