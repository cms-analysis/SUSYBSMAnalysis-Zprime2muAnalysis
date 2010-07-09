import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import Zprime2muAnalysisCommon

HistosFromPAT = cms.EDAnalyzer('Zprime2muHistosFromPAT',
                               Zprime2muAnalysisCommon,
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               useAllLeptons = cms.bool(True),
                               leptonsFromDileptons = cms.bool(False),
                               lowerMassWin = cms.double(0.5),
                               upperMassWin = cms.double(100.0),
                               )
