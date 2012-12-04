import FWCore.ParameterSet.Config as cms

# A filter for post-tuple filtering on the goodData results as stored
# in a TriggerResults object instead of filtering at tuple-making
# time.
import HLTrigger.HLTfilters.hltHighLevel_cfi
goodDataFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
goodDataFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
goodDataFilter.HLTPaths = ['goodDataAll'] # can set to just 'goodDataPrimaryVertexFilter', for example
goodDataFilter.andOr = False # = AND

from MuonPhotonMatch_cff import muonPhotonMatch
from OurSelectionDec2012_cff import allDimuons, dimuons, loose_cut

leptons = cms.EDProducer('Zprime2muLeptonProducer',
                         muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), #JMTBAD changeme after new PAT tuples
                         electron_src = cms.InputTag('cleanPatElectrons'),
                         muon_cuts = cms.string(loose_cut),
                         electron_cuts = cms.string('userInt("HEEPId") == 0'),
                         muon_track_for_momentum = cms.string('TunePNew'),
                         muon_photon_match_src = cms.InputTag('muonPhotonMatch'),
                         electron_muon_veto_dR = cms.double(-1),
                         trigger_match_max_dR = cms.double(0.2),
                         trigger_summary_src = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
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

def switch_to_old_selection(process):
    # Unfortunate names: what used to be New is now old
    import OurSelectionNew_cff
    process.leptons.muon_cuts = OurSelectionNew_cff.loose_cut
    process.allDimuons = OurSelectionNew_cff.allDimuons
    process.dimuons = OurSelectionNew_cff.dimuons

def switch_hlt_process_name(process, name):
    # JMTBAD better place for this fcn
    class warn_hlt_visitor(object):
        def enter(self, visitee):
            for attr in dir(visitee):
                if str(getattr(visitee, attr)) == 'HLT':
                    print 'warning: visitee %s, attribute %s has value HLT' % (visitee, attr)
        def leave(self, visitee):
            pass
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    for path_name, path in process.paths.iteritems(): # why does values() throw an exception?
        for label in ['TriggerResults', 'hltL1GtObjectMap', 'hltTriggerSummaryAOD']:
            old = cms.InputTag(label, '', 'HLT')
            new = cms.InputTag(label, '', name)
            massSearchReplaceAnyInputTag(path, old, new, verbose=False)
        path.visit(warn_hlt_visitor())
                                        
