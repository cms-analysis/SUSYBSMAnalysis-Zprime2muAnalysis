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

# A filter for post-tuple filtering on the goodData results as stored
# in a TriggerResults object instead of filtering at tuple-making
# time.
goodDataFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
goodDataFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
goodDataFilter.HLTPaths = ['goodDataAll'] # can set to just 'goodDataPrimaryVertexFilter', for example
goodDataFilter.andOr = False

from MuonPhotonMatch_cff import muonPhotonMatch
#from PATCandViewShallowCloneCombiner_cfi import allDimuons
from VBTFSelection_cff import vbtf_loose, allDimuons

leptons = cms.EDProducer('Zprime2muLeptonProducer',
                         muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'),
                         electron_src = cms.InputTag('cleanPatElectrons'),
                         muon_cuts = cms.string(vbtf_loose),
                         electron_cuts = cms.string('userInt("HEEPId") == 0'),
                         muon_track_for_momentum = cms.string('pmc'),
                         muon_photon_match_src = cms.InputTag('muonPhotonMatch')
                         )

dimuons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDimuons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(1),
                         back_to_back_cos_angle_min = cms.double(0.02),
                         vertex_chi2_max = cms.double(10),
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
