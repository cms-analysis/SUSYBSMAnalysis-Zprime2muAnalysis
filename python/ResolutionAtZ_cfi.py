import FWCore.ParameterSet.Config as cms

#from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction,hardInteraction_MiniAOD

ResolutionAtZ = cms.EDAnalyzer('ResolutionAtZ',
                                   #hardInteraction = hardInteraction,
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   )


