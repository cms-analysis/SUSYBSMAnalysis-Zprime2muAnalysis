import FWCore.ParameterSet.Config as cms

muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRDPtLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0))'),
    # Follow VBTF's matching criteria.
    maxDPtRel             = cms.double(1),
    maxDeltaR             = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)
)

# For the trigger match, currently HLT_Mu40_v5 is the lowest-pT
# unprescaled single muon path. Neither Spring11 nor Summer11 MC has
# such a path; emulate with Mu15, which is common to both. (This would
# not work if the MC paths ever have prescales...)
#
# Obviously, these values must be changed or at least checked when new
# trigger menus are implemented and used, in either MC or data.
trigger_pt_threshold = 40
offline_pt_threshold = 45
mc_trigger_paths = ['HLT_Mu15_v1', 'HLT_Mu15_v2']
trigger_paths = ['HLT_Mu40_v1', 'HLT_Mu40_v2', 'HLT_Mu40_v3', 'HLT_Mu40_v4', 'HLT_Mu40_v5']
old_trigger_paths = ['HLT_Mu30_v1', 'HLT_Mu30_v2'] # In runs 160329-163869, there was no HLT_Mu40 in the trigger menu. The next run with data, 165071, uses a trigger menu with HLT_Mu30_v3 and HLT_Mu40_v1. CheckPrescales downstream should behave correctly because of this configuration.

trigger_match =  ['(!triggerObjectMatchesByPath("%s",1,0).empty() && triggerObjectMatchesByPath("%s",1,0).at(0).pt() > %i)' % (n,n, trigger_pt_threshold) for n in mc_trigger_paths + old_trigger_paths]
trigger_match += ['!triggerObjectMatchesByPath("%s",1,0).empty()' % n for n in trigger_paths]
trigger_match = '(' + ' || '.join(trigger_match) + ')'

prescaled_trigger_pt_threshold = 15
prescaled_offline_pt_threshold = 20
prescaled_trigger_paths = ['HLT_Mu15_v1', 'HLT_Mu15_v2', 'HLT_Mu15_v3', 'HLT_Mu15_v4', 'HLT_Mu15_v5', 'HLT_Mu15_v6', 'HLT_Mu15_v7', 'HLT_Mu15_v8']
overall_prescale = 1080

prescaled_trigger_match = ' || '.join('!triggerObjectMatchesByPath("%s",1,0).empty()' % n for n in prescaled_trigger_paths)
