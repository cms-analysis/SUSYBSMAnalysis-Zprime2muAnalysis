import FWCore.ParameterSet.Config as cms

Zprime2muAnalysisCommon = cms.PSet(
    maxDileptons      = cms.uint32(1),
    doingElectrons    = cms.bool(False),
    useGen            = cms.bool(False),
    useSim            = cms.bool(False),
    useTrigger        = cms.bool(True),
    useRaw            = cms.bool(False),
    useReco           = cms.bool(True),
    dateHistograms    = cms.untracked.bool(True),
)

# By putting it in the analysis path, this module can be used to
# filter out events that do not pass our trigger selection, which are
# currently the OR of the highest-pT single muon trigger and double
# muon triggers. It does not go in Zprime2muAnalysisSequence by
# default; users must specifically include it.
import HLTrigger.HLTfilters.hltHighLevel_cfi
hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltFilter.HLTPaths = ['HLT_Mu9', 'HLT_DoubleMu3']
hltFilter.andOr = True # == OR

leptons = cms.EDProducer('Zprime2muLeptonProducer',
                         muon_src = cms.InputTag('cleanPatMuons'),
                         electron_src = cms.InputTag('cleanPatElectrons'),
                         muon_cuts = cms.string('pt > 20. && isolationR03.sumPt < 10'),
                         electron_cuts = cms.string(''),
                         muon_track_for_momentum = cms.string('pmc'),
                         muon_photon_match_src = cms.InputTag('muonPhotonMatch')
                         )

# JMTBAD muon-photon matching is done here using the default muon
# momentum and not whichever refit momentum will be eventually used in
# the analysis. For most muons this shouldn't make much difference as
# the matching is in eta-phi, but if the default is really bad and the
# selected refit manages to recover, the photon match may be
# wrong. The photon brem recovery needs to be studied more anyway.
muonPhotonMatch = cms.EDProducer('TrivialDeltaRViewMatcher',
                                 src     = cms.InputTag('cleanPatMuons'),
                                 matched = cms.InputTag('cleanPatPhotons'),
                                 distMin = cms.double(0.1)
                                 )

allDimuons = cms.EDProducer('PATCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string('')
                            )

dimuons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDimuons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(1)
                         )

#cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == -2'), # e.g. to select only mu+e- when the CandCombiner above made both mu+e- and mu-e+

Zprime2muAnalysisSequence = cms.Sequence(muonPhotonMatch * leptons * allDimuons * dimuons)

def rec_levels(process, new_track_types):
    process.leptons.muon_tracks_for_momentum = cms.vstring(*new_track_types)
    process.Zprime2muAnalysisSequence = cms.Sequence(process.muonPhotonMatch * process.leptons)

    for t in new_track_types:
        ad = process.allDimuons.clone()
        label = 'leptons:%s' % t
        ad.decay = '%s@+ %s@-' % (label, label)
        setattr(process, 'allDimuons' + t, ad)

        d = process.dimuons.clone()
        d.src = 'allDimuons' + t
        setattr(process, 'dimuons' + t, d)

        process.Zprime2muAnalysisSequence *= ad
        process.Zprime2muAnalysisSequence *= d
