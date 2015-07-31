import FWCore.ParameterSet.Config as cms

HistosFromPAT_miniAOD = cms.EDAnalyzer('Zprime2muHistosFromPAT_miniAOD',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               leptonsFromDileptons = cms.bool(False),
                               beamspot_src = cms.InputTag('offlineBeamSpot'),
                               vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                               use_bs_and_pv = cms.bool(True),
                               )

HistosFromPAT = cms.EDAnalyzer('Zprime2muHistosFromPAT',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               leptonsFromDileptons = cms.bool(False),
                               beamspot_src = cms.InputTag('offlineBeamSpot'),
                               vertex_src = cms.InputTag('offlinePrimaryVertices'),
                               use_bs_and_pv = cms.bool(True),
                               )
