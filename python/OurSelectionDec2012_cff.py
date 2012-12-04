import FWCore.ParameterSet.Config as cms

# The high-pT muon selection updated at the end of 2012.
# Designed to be used ONLY with the updated TuneP momentum assignment.
#
# The starting point is the Muon POG tight muon selection, documented
# at
#
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonRecoPerformance2010
#
# We add cuts on isolation, on the 3D angle between muons to suppress
# cosmics, on the vertex chi2, and on dpT/pT to reject grossly
# mismeasured tracks (the latter three being implemented in
# Zprime2muCompositeCandidatePicker, since they can't be done by the
# StringCutParser). The dimuons must be opposite-sign, and in events
# with more than one dimuon passing all the cuts, we keep the highest
# mass one.
#
# So, both muons must pass this selection:
#
# - muon must be a global muon and a tracker muon (isGlobalMuon && isTrackerMuon)
# - pT > offline_pt_threshold
# - |dxy wrt beamspot| < 0.2 cm (abs(dB) < 0.2)
# - relative tracker isolation less than 10% (isolationR03.sumPt / innerTrack.pt < 0.10)
# - number of tracker layers with hits > 5 (globalTrack.hitPattern.trackerLayersWithMeasurement > 5)
# - at least one pixel hit (globalTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit; this implies trackerMuon (numberOfMatchedStations > 1)
# - dpT/pT < 0.3
#
# Then at least one muon must be trigger-matched to the single muon
# HLT path. (The single muon HLT path used will change as the trigger
# menu evolves with luminosity, so the details are kept in a single
# file as they get used multiple places.)
#
# So we have a loose-tight combiner that requires both muons to pass
# the above cuts ("loose" is then a misnomer), and at least one must
# pass the trigger match requirement (the only "tight" cut).

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, offline_pt_threshold

loose_cut = 'isGlobalMuon && ' \
            'isTrackerMuon && ' \
            'pt > %s && ' \
            'abs(dB) < 0.2 && ' \
            'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
            'globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatchedStations > 1'

loose_cut = loose_cut % offline_pt_threshold

tight_cut = trigger_match

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
                         do_remove_overlap = cms.bool(True),
                         back_to_back_cos_angle_min = cms.double(-0.9998), # this corresponds to the angle (pi - 0.02) rad = 178.9 deg
                         vertex_chi2_max = cms.double(10),
                         dpt_over_pt_max = cms.double(0.3)
                         )
