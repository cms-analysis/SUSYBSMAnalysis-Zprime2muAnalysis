import FWCore.ParameterSet.Config as cms

muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRDPtLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('type("TriggerMuon") && (path("HLT_Mu9*") || path("HLT_Mu15*") || path("HLT_Mu24*") || path("HLT_Mu30*") || path("HLT_Mu40*"))'),
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

trigger_match =  ['(!triggerObjectMatchesByPath("%s").empty() && triggerObjectMatchesByPath("%s").at(0).pt() > %i)' % (n,n, trigger_pt_threshold) for n in mc_trigger_paths]
trigger_match += ['!triggerObjectMatchesByPath("%s").empty()' % n for n in trigger_paths]
trigger_match = '(' + ' || '.join(trigger_match) + ')'

prescaled_trigger_pt_threshold = 15
prescaled_offline_pt_threshold = 20
prescaled_trigger_paths = ['HLT_Mu15_v1', 'HLT_Mu15_v2', 'HLT_Mu15_v3', 'HLT_Mu15_v4', 'HLT_Mu15_v5', 'HLT_Mu15_v6', 'HLT_Mu15_v7', 'HLT_Mu15_v8']
overall_prescale = 720

prescaled_trigger_match = ' || '.join('!triggerObjectMatchesByPath("%s").empty()' % n for n in prescaled_trigger_paths)
