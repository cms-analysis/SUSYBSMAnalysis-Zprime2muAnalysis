import FWCore.ParameterSet.Config as cms

# This is our old selection, as used in the 2010 limits paper.
#
# The starting point is the (old) VBTF selection, which was documented
# at (dig in the history)
#
# https://twiki.cern.ch/twiki/bin/view/CMS/VbtfZMuMuBaselineSelection
#
# We remove the cuts on muon pseudo-rapidity, change the isolation
# cut, and add cuts on the 3D angle between muons and the common
# vertex chi2 probability (the latter two being implemented in
# Zprime2muCompositeCandidatePicker, since they can't be done by the
# StringCutParser).
#
# So, both muons must pass this selection:
#
# - muon must be a global muon (isGlobalMuon)
# - cocktail pT > offline_pt_threshold
# - number of tracker hits >= 10 (innerTrack.hitPattern.numberOfValidTrackerHits >= 10)
# - relative tracker isolation less than 10% (isolationR03.sumPt / innerTrack.pt < 0.10)
#
# Then one muon must pass the above plus a tighter set of cuts:
#
# - |dxy wrt beamspot| < 0.2 cm (abs(dB) < 0.2)
# - muon global track chi2/ndof < 10 (globalTrack.normalizedChi2 < 10)
# - at least one pixel hit (innerTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit (globalTrack.hitPattern.muonStationsWithValidHits >= 2)
# - must be a tracker muon (isTrackerMuon)
# - trigger matching to the single muon HLT path
#
# (The single muon HLT path used will change as the trigger menu
# evolves with luminosity; the details are kept in another file for
# multiple uses.)
#
# So we have a combiner that requires
# both muons to pass the loose cut, and at least one must pass the
# tight cut.

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, offline_pt_threshold

loose_cut = 'isGlobalMuon && ' \
            'pt > %s && ' \
            'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
            'innerTrack.hitPattern.numberOfValidTrackerHits >= 10'

loose_cut = loose_cut % offline_pt_threshold

tight_cut = 'abs(dB) < 0.2 && ' \
            'globalTrack.normalizedChi2 < 10 && ' \
            'innerTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
            'globalTrack.hitPattern.muonStationsWithValidHits >= 2 && ' \
            'isTrackerMuon && ' + trigger_match

allDimuons = cms.EDProducer('Zprime2muCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string(''),
                            loose_cut = cms.string(loose_cut),
                            tight_cut = cms.string(tight_cut)
                            )

dimuons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDimuons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(1),
                         sort_by_pt = cms.bool(False),
                         do_remove_overlap = cms.bool(True),
                         back_to_back_cos_angle_min = cms.double(-0.9998), # this corresponds to the angle (pi - 0.02) rad = 178.9 deg
                         vertex_chi2_max = cms.double(10),
                         )
