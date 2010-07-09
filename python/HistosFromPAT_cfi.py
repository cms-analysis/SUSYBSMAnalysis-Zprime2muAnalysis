import FWCore.ParameterSet.Config as cms

HistosFromPAT = cms.EDAnalyzer('Zprime2muHistosFromPAT',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               useAllLeptons = cms.bool(True),
                               leptonsFromDileptons = cms.bool(False),
                               lowerMassWin = cms.double(0.5),
                               upperMassWin = cms.double(100.0),
                               )
