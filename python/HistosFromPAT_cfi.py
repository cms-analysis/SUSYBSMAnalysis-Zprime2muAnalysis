import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction,hardInteraction_MiniAOD

HistosFromPAT = cms.EDAnalyzer('Zprime2muHistosFromPAT',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               leptonsFromDileptons = cms.bool(False),
                               beamspot_src = cms.InputTag('offlineBeamSpot'),
                               vertex_src = cms.InputTag('offlinePrimaryVertices'),
                               use_bs_and_pv = cms.bool(True),
                               useMadgraphWeight = cms.bool(True),
                               usekFactor = cms.bool(False),
                               hardInteraction = hardInteraction,
                               )

HistosFromPAT_MiniAOD = cms.EDAnalyzer('Zprime2muHistosFromPAT',
                               lepton_src = cms.InputTag('leptons', 'muons'),
                               dilepton_src = cms.InputTag('dimuons'),
                               leptonsFromDileptons = cms.bool(False),
                               beamspot_src = cms.InputTag('offlineBeamSpot'),
                               vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                               use_bs_and_pv = cms.bool(True),
                               useMadgraphWeight = cms.bool(True),
                               usekFactor = cms.bool(False),
                               hardInteraction = hardInteraction_MiniAOD,
)
