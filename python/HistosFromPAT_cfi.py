import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction,hardInteraction_MiniAOD
from SUSYBSMAnalysis.Zprime2muAnalysis.LRWeightProducer_cfi import LRWeightProducer

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
                               pu_src = cms.InputTag('slimmedAddPileupInfo'),
                               use_bs_and_pv = cms.bool(True),
                               useMadgraphWeight = cms.bool(True),
                               usekFactor = cms.bool(False),
                               useTTBarWeight = cms.bool(False),
			       hardInteraction = hardInteraction_MiniAOD,
                               doElectrons = cms.bool(False),
                               pu_weights = cms.vstring(),
			       year = cms.int32(2017),
              		       lrWeightProducer = LRWeightProducer,
			       prefireWeights = cms.InputTag("prefiringweight:nonPrefiringProb"),
			       LHEInfo = cms.InputTag("source")	 
)
