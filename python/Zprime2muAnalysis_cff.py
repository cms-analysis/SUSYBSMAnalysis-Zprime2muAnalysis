import FWCore.ParameterSet.Config as cms

# A filter for post-tuple filtering on the goodData results as stored
# in a TriggerResults object instead of filtering at tuple-making
# time.
import HLTrigger.HLTfilters.hltHighLevel_cfi
goodDataFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
goodDataFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
goodDataFilter.HLTPaths = ["goodDataPrimaryVertexFilter"] # can set to just 'goodDataPrimaryVertexFilter', for example
#goodDataFilter.HLTPaths = ['goodDataMETFilter']
goodDataFilter.andOr = False # = AND

from SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff import primaryVertexMiniAOD, hltPhysicsDeclared, metFilters
goodDataFiltersMiniAOD = [primaryVertexMiniAOD]
## for full filtering, use:
#goodDataFiltersMiniAOD = [primaryVertexMiniAOD,hltPhysicsDeclared]
#goodDataFiltersMiniAOD += metFilters


from MuonPhotonMatch_cff import muonPhotonMatch, muonPhotonMatchMiniAOD
from OurSelection2016_cff import allDimuons, dimuons, loose_cut
#from OurSelectionDec2012_cff import allDimuons, dimuons, loose_cut

leptons = cms.EDProducer('Zprime2muLeptonProducer',
                         muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), #JMTBAD changeme after new PAT tuples
                         electron_src = cms.InputTag('cleanPatElectrons'),
                         muon_srcSecond = cms.InputTag('cleanPatMuonsTriggerMatch'), #JMTBAD changeme after new PAT tuples
                         electron_srcSecond = cms.InputTag('cleanPatElectrons'),
                         muon_cuts = cms.string(loose_cut),
                         electron_cuts = cms.string(''),
                         muon_track_for_momentum = cms.string('TunePNew'),
                         muon_track_for_momentum_CSC = cms.string('Inner'),
                         muon_photon_match_src = cms.InputTag('muonPhotonMatch'),
                         electron_muon_veto_dR = cms.double(-1),
                         trigger_match_max_dR = cms.double(0.2),
                         trigger_summary_src = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
                         )
leptonsMini = cms.EDProducer('Zprime2muLeptonProducer_miniAOD',
                              muon_src = cms.InputTag('slimmedMuons'), #JMTBAD changeme after new PAT tuples
                              electron_src = cms.InputTag('slimmedElectrons'),
                              electron_id = cms.InputTag('egmGsfElectronIDs:heepElectronID-HEEPV60'),
                              muon_srcSecond = cms.InputTag('slimmedMuons'), #JMTBAD changeme after new PAT tuples
                              muon_cuts = cms.string(loose_cut),
                              muon_track_for_momentum = cms.string('TunePNew'),
                              muon_track_for_momentum_CSC = cms.string('Inner'),
                              muon_photon_match_src = cms.InputTag('muonPhotonMatchMiniAOD'),
                              electron_muon_veto_dR = cms.double(-1),
                              trigger_match_max_dR = cms.double(0.2),
                              trigger_summary = cms.InputTag('selectedPatTrigger'),
#                              bits = cms.InputTag("TriggerResults","","HLT2"),##mc reHLT
                              bits = cms.InputTag("TriggerResults","","HLT"),#data
                              prescales = cms.InputTag("patTrigger"),
                              )

Zprime2muAnalysisSequence = cms.Sequence(muonPhotonMatch * leptons * allDimuons * dimuons)
Zprime2muAnalysisSequence_MiniAOD = cms.Sequence(muonPhotonMatchMiniAOD * leptonsMini * allDimuons * dimuons)


#####################################################################
############# E L E C T R O N  -  S E L E C T O R ###################
# defined a proces to call VID Ele Selector using miniAOD 
# the proces is called in test/DataMCSpectraComparison/histos.py
# and in test/NMinus1Effs/nminus1effs.py

# to use it you have to do from CMSSW_X_X_X/src/:
# git cms-merge-topic Sam-Harper:HEEPV70VID
# compile it, and pass to "my_id_modules" the version that you merged 

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
def electrons_miniAOD(process):
    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    # my_id_module is the 
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#####################################################################    

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
        #path.visit(warn_hlt_visitor())
 
def switch_reco_process_name(process, name):
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
        for label in ['TriggerResults']:
            old = cms.InputTag(label, '', 'RECO')
            new = cms.InputTag(label, '', name)
            massSearchReplaceAnyInputTag(path, old, new, verbose=False)
        #path.visit(warn_hlt_visitor())
                                        
