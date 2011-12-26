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

# For the trigger match, currently HLT_Mu40_eta2p1_v1 is the lowest-pT
# unprescaled single muon path, albeit with coverage in |eta| only up
# to 2.1. Neither Spring11 nor Summer11 MC has such a path; emulate
# with Mu15, which is common to both. (This would not work if the MC
# paths ever have prescales...)
#
# Obviously, these values must be changed or at least checked when new
# trigger menus are implemented and used, in either MC or data.

trigger_pt_threshold = 40
offline_pt_threshold = 45

mc_trigger_paths = ['HLT_Mu15_v1', 'HLT_Mu15_v2']
trigger_paths = ['HLT_Mu40_eta2p1_v%i' % i for i in (1,4,5)]
old_trigger_paths = ['HLT_Mu40_v%i' % i for i in (1,2,3,5)] + ['HLT_Mu30_v1', 'HLT_Mu30_v2']

trigger_match =  ['(!triggerObjectMatchesByPath("%(path)s",1,0).empty() && triggerObjectMatchesByPath("%(path)s",1,0).at(0).pt() > %(trigger_pt_threshold)i && abs(triggerObjectMatchesByPath("%(path)s",1,0).at(0).eta()) < 2.1)' % locals() for path in mc_trigger_paths + old_trigger_paths]
trigger_match += ['!triggerObjectMatchesByPath("%s",1,0).empty()' % n for n in trigger_paths]
trigger_match = '(' + ' || '.join(trigger_match) + ')'

prescaled_trigger_pt_threshold = 15
prescaled_offline_pt_threshold = 20
prescaled_trigger_paths = ['HLT_Mu15_v%i' % i for i in (2,3,4,5,6,8,9,12,13)]
overall_prescale = 2000

prescaled_trigger_match = ' || '.join('(!triggerObjectMatchesByPath("%(path)s",1,0).empty() && abs(triggerObjectMatchesByPath("%(path)s",1,0).at(0).eta()) < 2.1)' % locals() for path in prescaled_trigger_paths)
