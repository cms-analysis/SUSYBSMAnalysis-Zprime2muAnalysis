import FWCore.ParameterSet.Config as cms

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

from MuonPhotonMatch_cff import muonPhotonMatch

allDimuons = cms.EDProducer('PATCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string('')
                            )

dimuons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDimuons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(1)
                         )

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
