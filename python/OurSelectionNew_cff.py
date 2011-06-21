import FWCore.ParameterSet.Config as cms

# The starting point is the Muon POG tight muon selection, documented
# at
#
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonRecoPerformance2010
#
# We add cuts on isolation, on the 3D angle between muons to suppress
# cosmics, and on the vertex chi2 (the latter two being implemented in
# Zprime2muCompositeCandidatePicker, since they can't be done by the
# StringCutParser). The dimuons must be opposite-sign, and in events
# with more than one dimuon passing all the cuts, we keep the highest
# mass one.
#
# So, both muons must pass this selection:
#
# - muon must be a global muon and a tracker muon (isGlobalMuon && isTrackerMuon)
# - pT > 35 (pt > 35.)
# - |dxy wrt beamspot| < 0.2 cm (abs(dB) < 0.2)
# - muon global track chi2/ndof < 10 (globalTrack.normalizedChi2 < 10)
# - relative tracker isolation less than 10% (isolationR03.sumPt / innerTrack.pt < 0.10)
# - number of tracker hits > 10 (globalTrack.hitPattern.numberOfValidTrackerHits > 10)
# - at least one pixel hit (globalTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit; this implies trackerMuon (numberOfMatchedStations > 1)
#
# Then at least one muon must be trigger-matched to the single muon
# HLT path (e.g. !triggerObjectMatchesByPath("HLT_Mu30").empty()) (The
# single muon HLT path used will change as the trigger menu evolves
# with luminosity.)
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the above cuts ("loose" is then a misnomer), and
# at least one must pass the trigger match requirement (the only
# "tight" cut).

loose_cut = 'isGlobalMuon && ' \
            'isTrackerMuon && ' \
            'pt > 35. && ' \
            'abs(dB) < 0.2 && ' \
            'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatchedStations > 1'

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match

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
