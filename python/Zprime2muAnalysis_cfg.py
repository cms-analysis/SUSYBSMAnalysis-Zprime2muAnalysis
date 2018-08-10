import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:patTuneP.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service('TFileService', fileName=cms.string('zp2mu_histos.root'))

#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'
#process.GlobalTag.globaltag = '76X_dataRun2_v15'
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff')
