import FWCore.ParameterSet.Config as cms

LRWeightProducer = cms.PSet(src = cms.InputTag('prunedMCLeptons'),
                           doingElectrons = cms.bool(False),
                           doingLR = cms.bool(True),
                           calculate = cms.bool(False),
			   interference = cms.int32(1),
			   Lambda = cms.int32(16000),
                           )


