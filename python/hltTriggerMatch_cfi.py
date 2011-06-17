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

# For the trigger match, currently HLT_Mu30_v3 is the lowest-pT
# unprescaled single muon path. Spring11 MC does not have such a path;
# emulate with Mu15.
trigger_match = '(' \
                '(!triggerObjectMatchesByPath("HLT_Mu15_v1").empty() && triggerObjectMatchesByPath("HLT_Mu15_v1").at(0).pt() > 30) || ' \
                '(!triggerObjectMatchesByPath("HLT_Mu15_v2").empty() && triggerObjectMatchesByPath("HLT_Mu15_v2").at(0).pt() > 30) || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v1").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v2").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v3").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v4").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v5").empty()' \
                ')'
