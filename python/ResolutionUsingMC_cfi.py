import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction

ResolutionUsingMC = cms.EDAnalyzer('ResolutionUsingMC',
                                   hardInteraction = hardInteraction,
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   useAllLeptons = cms.bool(True),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   lowerMassWin = cms.double(0.5),
                                   upperMassWin = cms.double(100.0),
                                   )
