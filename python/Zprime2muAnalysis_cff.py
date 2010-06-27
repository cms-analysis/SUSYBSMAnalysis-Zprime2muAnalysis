import FWCore.ParameterSet.Config as cms

Zprime2muAnalysisCommon = cms.PSet(
    ################################################################
    # Analysis configuration:
    ################################################################
    verbosity         = cms.untracked.int32(0),
    maxDileptons      = cms.uint32(1),
    doingElectrons    = cms.bool(False),
    useGen            = cms.bool(False),
    useSim            = cms.bool(False),
    useTrigger        = cms.bool(True),
    useRaw            = cms.bool(False),
    useReco           = cms.bool(True),
    dateHistograms    = cms.untracked.bool(True),
    cutMask           = cms.uint32(0),
    
    ################################################################
    # Trigger info: InputTags and path names.
    ################################################################
    l1GtObjectMap = cms.InputTag('hltL1GtObjectMap'),
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # Trigger paths we use: the result is the OR of these. These
    # are the highest pT muon triggers in the 8E29 menu. The 1E31
    # menu has in addition L1_SingleMu10 which seeds HLT_Mu15.
    l1Paths = cms.vstring('L1_SingleMu7', 'L1_DoubleMu3'),
    hltPaths = cms.vstring('HLT_Mu9', 'HLT_DoubleMu3'),
)

muon_track_types = ['global', 'inner', 'outer', 'tpfms', 'picky', 'pmc', 'tmr', 'sigmaswitch']

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
