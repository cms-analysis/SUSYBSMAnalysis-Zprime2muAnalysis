import FWCore.ParameterSet.Config as cms

# The starting point is the VBTF selection, which is documented at
#
# https://twiki.cern.ch/twiki/bin/view/CMS/VbtfZMuMuBaselineSelection
#
# We remove the cuts on muon pseudo-rapidity, change the isolation cut
# to 10 GeV, and add cuts on the 3D angle between muons and the common
# vertex chi2 probability (the latter two being implemented in
# Zprime2muCompositeCandidatePicker, since they can't be done by the
# StringCutParser).
#
# So, both muons must pass this selection:
#
# - muon must be a global muon
# - pT > 20
# - number of tracker hits >= 10 (innerTrack.hitPattern.numberOfValidTrackerHits >= 10)
#
# Then one muon must pass the above plus a tighter set of cuts:
#
# - dxy < 0.2 cm (dB < 0.2)
# - muon global track chi2/ndof < 10 (globalTrack.normalizedChi2 < 10)
# - at least one pixel hit (innerTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit (globalTrack.hitPattern.muonStationsWithValidHits >= 2)
# - must be a tracker muon
# - trigger matching to the single muon HLT path (e.g. !triggerObjectMatchesByPath("HLT_Mu9").empty())
#
# (The single muon HLT path used will change as the trigger menu
# evolves with luminosity.)
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the loose cut, and at least one must pass the
# tight cut.

vbtf_loose = 'isGlobalMuon && ' \
             'pt > 20. && ' \
             'isolationR03.sumPt < 10 &&' \
             'innerTrack.hitPattern.numberOfValidTrackerHits >= 10'

# For the trigger match, currently HLT_Mu15_v1 is the lowest-pT
# unprescaled single muon path. In runs <= 147119, HLT_Mu15_v1 did not
# exist. Emulate it by using HLT_Mu9 (unprescaled in those runs) and a
# pT cut.
vbtf_trigger_match = '(' \
                     '(!triggerObjectMatchesByPath("HLT_Mu9").empty() && triggerObjectMatchesByPath("HLT_Mu9").at(0).pt() > 15) || ' \
                     '!triggerObjectMatchesByPath("HLT_Mu15_v1").empty()' \
                     ')'

vbtf_tight = 'dB < 0.2 && ' \
             'globalTrack.normalizedChi2 < 10 && ' \
             'innerTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
             'globalTrack.hitPattern.muonStationsWithValidHits >= 2 && ' \
             'isTrackerMuon && ' + vbtf_trigger_match

allDimuons = cms.EDProducer('LooseTightCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string(''),
                            loose_cut = cms.string(vbtf_loose),
                            tight_cut = cms.string(vbtf_tight)
                            )
