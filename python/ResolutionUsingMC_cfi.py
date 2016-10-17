import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction,hardInteraction_MiniAOD

ResolutionUsingMC = cms.EDAnalyzer('ResolutionUsingMC',
                                   hardInteraction = hardInteraction,
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   )

ResolutionUsingMC_MiniAOD = cms.EDAnalyzer('ResolutionUsingMC',
                                   hardInteraction = hardInteraction_MiniAOD,
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   )
