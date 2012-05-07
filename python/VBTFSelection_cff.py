import FWCore.ParameterSet.Config as cms
#
# This attempts to implement exactly the VBTF selection, for which the
# twiki is out of date. The most up-to-date info comes from AN-10-395 and
#
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/ElectroWeakAnalysis/ZMuMu/python/ZMuMuGolden_cfi.py?revision=1.7&view=markup
#
# (Since the single muon trigger threshold is now above the old one,
# e.g. pT > 40 GeV, here the offline pT cut is raised.)
#
# Their golden dimuons are formed from opposite-charged
# pairs in which both muons must pass this selection:
#
# - muon must be a global muon (isGlobalMuon)
# - muon must be a tracker muon (isTrackerMuon)
# - pT > offline_pt_threshold
# - |eta| < 2.1 (abs(innerTrack.eta) < 2.1)
# - |dxy wrt beamspot| < 0.2 cm (abs(dB) < 0.2)
# - absolute tracker-only isolation < 3 GeV (isolationR03.sumPt < 3)
# - number of tracker hits > 10 (innerTrack.hitPattern.numberOfValidTrackerHits > 10)
# - at least one pixel hit (innerTrack.hitPattern.numberOfValidPixelHits > 0)
# - at least two muon segments matched (numberOfMatches > 1)
#
# The ZMuMuGolden_cfi.py linked above uses abs(globalTrack.dxy) < 1.0
# for some reason. Also it has a cut on the global-track chi2, which
# AN-10-395 says isn't used. Follow
#
# Then at least one muon must be trigger-matched to the single muon
# HLT path. (The single muon HLT path used will change as the trigger
# menu evolves with luminosity, so the details are kept in a single
# file as they get used multiple places.)
#
# So we have a loose-tight combiner that requires
# both muons to pass the above cuts ("loose" is then a misnomer), and
# at least one must pass the trigger match requirement (the only
# "tight" cut).

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, offline_pt_threshold

loose_cut = 'isGlobalMuon && ' \
            'isTrackerMuon && ' \
            'innerTrack.pt > %s && ' \
            'abs(innerTrack.eta) < 2.1 && ' \
            'isolationR03.sumPt < 3 && ' \
            'abs(dB) < 0.2 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits > 0 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatches > 1'

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
                         max_candidates = cms.uint32(100),
                         do_remove_overlap = cms.bool(False),
                         )
