import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("METFilterMiniAOD",
                               src = cms.InputTag("TriggerResults","","PAT"),	# MC
#                                src = cms.InputTag("TriggerResults","","RECO"), # data
                               flag = cms.string("Flag_globalTightHalo2016Filter"),
                               applyfilter = cms.untracked.bool(True),
                               debugOn = cms.untracked.bool(False),
                               
                               )
