import FWCore.ParameterSet.Config as cms

# https://twiki.cern.ch/twiki/bin/view/CMS/VbtfZMuMuBaselineSelection
#
# Both muons must pass this selection:
#
# - muon must be a global muon
# - pT > 20
# - |eta| < 2.4
# - tracker isolation < 3 GeV (isolationR03.sumPt < 3)
# - number of tracker hits >= 10 (innerTrack.hitPattern.numberOfValidTrackerHits >= 10)
#
# Then one muon must pass the above plus a tighter set of cuts:
#
# - dxy < 0.2 cm (dB < 0.2)
# - muon global track chi2/ndof < 10 (globalTrack.normalizedChi2 < 10)
# - at least one pixel hit (innerTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit (globalTrack.hitPattern.muonStationsWithValidHits >= 2)
# - must be a tracker muon
# - |eta| < 2.1
# - trigger matching to the single muon HLT path (e.g. !triggerObjectMatchesByPath("HLT_Mu9").empty())
# - the L3 muon so matched has |eta| < 2.1 (e.g. abs(triggerObjectMatchesByPath("HLT_Mu9").at(0).eta()) < 2.1)
#
# (The single muon HLT path used will change as the trigger menu
# evolves with luminosity.)
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the loose cut, and at least one must pass the
# tight cut.

vbtf_loose = 'isGlobalMuon && ' \
             'pt > 20. && ' \
             'abs(eta) < 2.4 && ' \
             'isolationR03.sumPt < 3 && ' \
             'innerTrack.hitPattern.numberOfValidTrackerHits >= 10'

# For the trigger match, currently HLT_Mu11 is the lowest-pT
# unprescaled single muon path. In runs <= 147119, HLT_Mu11 did not
# exist. Emulate it by using HLT_Mu9 (unprescaled in those runs) and a
# pT cut.
vbtf_trigger_match = '(' \
                     '(!triggerObjectMatchesByPath("HLT_Mu9").empty() && abs(triggerObjectMatchesByPath("HLT_Mu9").at(0).eta()) < 2.1 && triggerObjectMatchesByPath("HLT_Mu9").at(0).pt() > 11) || ' \
                     '(!triggerObjectMatchesByPath("HLT_Mu11").empty() && abs(triggerObjectMatchesByPath("HLT_Mu11").at(0).eta()) < 2.1)' \
                     ')'

vbtf_tight = 'dB < 0.2 && ' \
             'globalTrack.normalizedChi2 < 10 && ' \
             'innerTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
             'globalTrack.hitPattern.muonStationsWithValidHits >= 2 && ' \
             'isTrackerMuon && ' \
             'abs(eta) < 2.1 && ' + vbtf_trigger_match

allDimuons = cms.EDProducer('LooseTightCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string(''),
                            loose_cut = cms.string(vbtf_loose),
                            tight_cut = cms.string(vbtf_tight)
                            )
