import FWCore.ParameterSet.Config as cms

CheckPrescale = cms.EDAnalyzer('CheckPrescale',
                               hlt_src = cms.InputTag('TriggerResults', '', 'HLT'),
                               trigger_paths = cms.vstring(),
                               dump_prescales = cms.untracked.bool(False),
                               )
