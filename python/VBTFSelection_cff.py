import FWCore.ParameterSet.Config as cms

# This attempts to implement exactly the VBTF selection, for which the
# twiki is out of date. The most up-to-date info comes from AN-10-395 and
#
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/ElectroWeakAnalysis/ZMuMu/python/ZMuMuGolden_cfi.py?revision=1.7&view=markup
#
# (Since the single muon trigger threshold is now pT > 30 GeV, here
# the offline pT cut is raised to 35 from 20 GeV.)
#
# Their golden dimuons are formed from opposite-charged
# pairs in which both muons must pass this selection:
#
# - muon must be a global muon (isGlobalMuon)
# - muon must be a tracker muon (isTrackerMuon)
# - pT > 35 GeV (innerTrack.pt > 35)
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
# At least one muon must be matched to a trigger object firing the
# single muon HLT path
# (e.g. !triggerObjectMatchesByPath("HLT_Mu9").empty()) The single
# muon HLT path used will change as the trigger menu evolves with
# luminosity.
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the above cuts ("loose" is then a misnomer), and
# at least one must pass the trigger match requirement (the only
# "tight" cut).

loose_cut = 'isGlobalMuon && ' \
            'isTrackerMuon && ' \
            'innerTrack.pt > 35. && ' \
            'abs(innerTrack.eta) < 2.1 && ' \
            'isolationR03.sumPt < 3 && ' \
            'abs(dB) < 0.2 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits > 0 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatches > 1'

# For the trigger match, currently HLT_Mu30_v3 is the lowest-pT
# unprescaled single muon path. Spring11 MC does not have such a path;
# emulate with Mu15.
trigger_match = '(' \
                '(!triggerObjectMatchesByPath("HLT_Mu15_v1").empty() && triggerObjectMatchesByPath("HLT_Mu15_v1").at(0).pt() > 30) || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v1").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v2").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu30_v3").empty()' \
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
                         max_candidates = cms.uint32(100),
                         do_remove_overlap = cms.bool(False),
                         )
