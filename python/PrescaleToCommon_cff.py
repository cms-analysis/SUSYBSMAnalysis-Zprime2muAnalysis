import FWCore.ParameterSet.Config as cms

RandomNumberGeneratorService = cms.Service('RandomNumberGeneratorService',
                                           PrescaleToCommon = cms.PSet(initialSeed = cms.untracked.uint32(1219))
                                           )

PrescaleToCommon = cms.EDFilter('PrescaleToCommon',
                                hlt_src = cms.InputTag('TriggerResults','','HLT'),
                                trigger_paths = cms.vstring(),
                                overall_prescale = cms.int32(1),
                                assume_simulation_has_prescale_1 = cms.bool(True) # Current PAT tuples of MC samples don't have both L1 branches :-(
                                )
