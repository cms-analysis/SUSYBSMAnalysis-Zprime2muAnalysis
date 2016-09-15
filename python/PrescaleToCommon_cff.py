import FWCore.ParameterSet.Config as cms

RandomNumberGeneratorService = cms.Service('RandomNumberGeneratorService',
                                           PrescaleToCommon = cms.PSet(initialSeed = cms.untracked.uint32(1219)),
                                           PrescaleToCommonMiniAOD = cms.PSet(initialSeed = cms.untracked.uint32(1219))
                                           )

PrescaleToCommon = cms.EDFilter('PrescaleToCommon',
                                hlt_src = cms.InputTag('TriggerResults','','HLT2'),
                                TriggerResults_src = cms.InputTag('TriggerResults', '', 'HLT2'),
                                trigger_paths = cms.vstring(),
                                overall_prescale = cms.int32(1),
                                assume_simulation_has_prescale_1 = cms.bool(True) # Current PAT tuples of MC samples don't have both L1 branches :-(
                                )

PrescaleToCommonMiniAOD = cms.EDFilter('PrescaleToCommon_miniAOD',
                                hlt_src = cms.InputTag('TriggerResults','','HLT2'),
                                TriggerResults_src = cms.InputTag('TriggerResults', '', 'HLT2'),
				Prescale_src = cms.InputTag('patTrigger','','RECO'),
				L1Prescale_max_src = cms.InputTag('patTrigger','l1max','RECO'),
				L1Prescale_min_src = cms.InputTag('patTrigger','l1min','RECO'),
                                trigger_paths = cms.vstring(),
                                overall_prescale = cms.int32(1),
                                assume_simulation_has_prescale_1 = cms.bool(True) # Current PAT tuples of MC samples don't have both L1 branches :-(                                
			        
)

