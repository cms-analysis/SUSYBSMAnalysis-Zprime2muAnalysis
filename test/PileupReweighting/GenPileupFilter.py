import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('/store/mc/Summer11/DYToMuMu_M-1000_TuneZ2_7TeV-pythia6-tauola/GEN-SIM-RECO/PU_S3_START42_V11-v2/0000/E8CBB8E9-8988-E011-A764-1CC1DE1CDF2A.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service('TFileService', fileName=cms.string('zp2mu_histos.root'))

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.GenPileupFilter_cfi')
process.p = cms.Path(process.GenPileupFilter)
