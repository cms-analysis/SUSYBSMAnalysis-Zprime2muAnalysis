import FWCore.ParameterSet.Config as cms

# The starting point is the (old) VBTF selection, which is documented
# at
#
# https://twiki.cern.ch/twiki/bin/view/CMS/VbtfZMuMuBaselineSelection
#
# We remove the cuts on muon pseudo-rapidity, change the isolation
# cut, and add cuts on the 3D angle between muons and the common
# vertex chi2 probability (the latter two being implemented in
# Zprime2muCompositeCandidatePicker, since they can't be done by the
# StringCutParser).
#
# Both muons must pass this selection:
#
# - muon must be a global muon and a tracker muon (isGlobalMuon && isTrackerMuon)
# - pT > 20 (innerTrack.pt > 20.)
# - |dxy wrt beamspot| < 0.2 cm (abs(dB) < 0.2)
# - relative tracker isolation less than 10% (isolationR03.sumPt / innerTrack.pt < 0.10)
# - number of tracker hits > 10 (globalTrack.hitPattern.numberOfValidTrackerHits > 10)
# - at least one pixel hit (globalTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit (globalTrack.hitPattern.muonStationsWithValidHits >= 2)
#
# Then at least one muon must be trigger-matched to the single muon
# HLT path (e.g. !triggerObjectMatchesByPath("HLT_Mu9").empty()) (The
# single muon HLT path used will change as the trigger menu evolves
# with luminosity.)
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the above cuts ("loose" is then a misnomer), and
# at least one must pass the trigger match requirement (the only
# "tight" cut).

loose_cut = 'isGlobalMuon && ' \
            'isTrackerMuon && ' \
            'innerTrack.pt > 20. && ' \
            'abs(dB) < 0.2 && ' \
            'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'globalTrack.hitPattern.muonStationsWithValidHits >= 2'

# For the trigger match, currently HLT_Mu15_v1 is the lowest-pT
# unprescaled single muon path. In runs <= 147119, HLT_Mu15_v1 did not
# exist. Emulate it by using HLT_Mu9 (unprescaled in those runs) and a
# pT cut.
trigger_match = '(' \
                '(!triggerObjectMatchesByPath("HLT_Mu9").empty() && triggerObjectMatchesByPath("HLT_Mu9").at(0).pt() > 15) || ' \
                '!triggerObjectMatchesByPath("HLT_Mu15_v1").empty()' \
                ')'

tight_cut = trigger_match

allDimuons = cms.EDProducer('LooseTightCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string(''),
                            loose_cut = cms.string(loose_cut),
                            tight_cut = cms.string(tight_cut)
                            )

dimuons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDimuons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(1),
                         do_remove_overlap = cms.bool(True),
                         back_to_back_cos_angle_min = cms.double(-0.9998), # this corresponds to the angle (pi - 0.02) rad = 178.9 deg
                         vertex_chi2_max = cms.double(10),
                         )
