import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import Zprime2muAnalysisCommon

ResolutionUsingMC = cms.EDAnalyzer('ResolutionUsingMC',
                                   Zprime2muAnalysisCommon.clone(useGen = True),
                                   lepton_src = cms.InputTag('leptons', 'muons'),
                                   dilepton_src = cms.InputTag('dimuons'),
                                   triggeredEventsOnly = cms.bool(True),
                                   useAllLeptons = cms.bool(True),
                                   leptonsFromDileptons = cms.bool(False),
                                   doQoverP = cms.bool(False),
                                   lowerMassWin = cms.double(0.5),
                                   upperMassWin = cms.double(100.0),
                                   )
