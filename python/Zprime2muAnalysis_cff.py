import FWCore.ParameterSet.Config as cms

# By putting it in the analysis path, this module can be used to
# filter out whole events that do not pass our trigger selection,
# which is currently the highest-pT unprescaled single muon trigger in
# the 2E32 menu, HLT_Mu15_v1. When using one of the two selections
# VBTFSelection and OurSelection, one muon is required to match to a
# trigger object passing the single muon trigger anyway, so using this
# is redundant in some cases. It does not go in
# Zprime2muAnalysisSequence by default; users must specifically
# include it.
import HLTrigger.HLTfilters.hltHighLevel_cfi
hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltFilter.HLTPaths = ['HLT_Mu15_v1']

# A filter for post-tuple filtering on the goodData results as stored
# in a TriggerResults object instead of filtering at tuple-making
# time.
goodDataFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
goodDataFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
goodDataFilter.HLTPaths = ['goodDataAll'] # can set to just 'goodDataPrimaryVertexFilter', for example
goodDataFilter.andOr = False # = AND

from MuonPhotonMatch_cff import muonPhotonMatch
#from PATCandViewShallowCloneCombiner_cfi import allDimuons
#from VBTFSelection_cff import allDimuons, dimuons, loose_cut
from OurSelection_cff import allDimuons, dimuons, loose_cut

leptons = cms.EDProducer('Zprime2muLeptonProducer',
                         muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'),
                         electron_src = cms.InputTag('cleanPatElectrons'),
                         muon_cuts = cms.string(loose_cut),
                         electron_cuts = cms.string('userInt("HEEPId") == 0'),
                         muon_track_for_momentum = cms.string('pmc'),
                         muon_photon_match_src = cms.InputTag('muonPhotonMatch')
                         )

Zprime2muAnalysisSequence = cms.Sequence(muonPhotonMatch * leptons * allDimuons * dimuons)

def rec_levels(process, new_track_types):
    process.leptons.muon_tracks_for_momentum = cms.vstring(*new_track_types)
    process.Zprime2muAnalysisSequence = cms.Sequence(process.muonPhotonMatch * process.leptons)
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.muonPhotonMatch * process.leptons * process.allDimuons * process.dimuons)

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

def rec_level_module(process, module, name, tracks):
    p = []
    for t in tracks:
        h = module.clone()
        if hasattr(h, 'lepton_src'):
            h.lepton_src = cms.InputTag('leptons', t)
        if hasattr(h, 'dilepton_src'):
            h.dilepton_src = cms.InputTag('dimuons' + t)
        setattr(process, name + t, h)
        p.append(h)
    return reduce(lambda x,y: x*y, p)
