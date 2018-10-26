import FWCore.ParameterSet.Config as cms

loose_cut = ''

tight_cut = ''

allDielectrons = cms.EDProducer('Zprime2muCombiner',
                            decay = cms.string('leptons:electrons@+ leptons:electrons@-'),
                            cut = cms.string(''),
                            loose_cut = cms.string(loose_cut),
                            tight_cut = cms.string(tight_cut)
                            )

dielectrons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDielectrons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(1),
                         sort_by_pt = cms.bool(True),
                         prefer_Z = cms.bool(True),
                         do_remove_overlap = cms.bool(True),
                         back_to_back_cos_angle_min = cms.double(-2), # this corresponds to the angle (pi - 0.02) rad = 178.9 deg
#                         vertex_chi2_max = cms.double(10),
                         vertex_chi2_max = cms.double(999999999999),
                         dpt_over_pt_max = cms.double(999999999999)
                         )
dielectronHLT = cms.EDFilter("TriggerResultsFilter",
    	triggerConditions = cms.vstring("HLT_DoubleEle33_CaloIdL_MW_v*"),
	hltResults = cms.InputTag("TriggerResults","","HLT"),
	l1tResults = cms.InputTag("gtStage2Digis")
	)
