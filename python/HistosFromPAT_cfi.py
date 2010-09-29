import FWCore.ParameterSet.Config as cms

HistosFromPAT = cms.EDAnalyzer('Zprime2muHistosFromPAT',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               leptonsFromDileptons = cms.bool(False),
                               )
