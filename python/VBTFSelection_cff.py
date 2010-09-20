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
# - trigger matching to HLT_Mu9
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the loose cut, and at least one must pass the
# tight cut.
#
# To use this instead of the default allDimuons, in your top level
# analysis cfg you can include this cff, e.g.:
#
# process.load('SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff')
# process.load('SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff')
# process.allDimuons = process.allDimuonsVBTF

vbtf_loose = 'isGlobalMuon && pt > 20. && abs(eta) < 2.4 && isolationR03.sumPt < 3 && innerTrack.hitPattern.numberOfValidTrackerHits >= 10'
vbtf_tight = 'dB < 0.2 && globalTrack.normalizedChi2 < 10 && innerTrack.hitPattern.numberOfValidPixelHits >= 1 && globalTrack.hitPattern.muonStationsWithValidHits >= 2 && isTrackerMuon && abs(eta) < 2.1 && !triggerObjectMatchesByPath("HLT_Mu9").empty()'

allDimuonsVBTF = cms.EDProducer('LooseTightCandViewShallowCloneCombiner',
                                decay = cms.string('leptons:muons@+ leptons:muons@-'),
                                cut = cms.string(''),
                                loose_cut = cms.string(vbtf_loose),
                                tight_cut = cms.string(vbtf_tight)
                                )

