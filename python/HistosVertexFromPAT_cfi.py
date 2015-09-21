import FWCore.ParameterSet.Config as cms

HistosVertexFromPAT = cms.EDAnalyzer('Zprime2muHistosVertexFromPAT',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               leptonsFromDileptons = cms.bool(False),
                               beamspot_src = cms.InputTag('offlineBeamSpot'),
                               vertex_src = cms.InputTag('offlinePrimaryVertices'),
                               use_bs_and_pv = cms.bool(True),
                               )
