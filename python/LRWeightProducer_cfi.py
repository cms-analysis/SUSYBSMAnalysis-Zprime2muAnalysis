import FWCore.ParameterSet.Config as cms

LRWeightProducer = cms.PSet(src = cms.InputTag('prunedMCLeptons'),
                           doingElectrons = cms.bool(False),
                           doingLR = cms.bool(True),
                           calculate = cms.bool(False),
			   Lambda = cms.int32(16000),
                           )


