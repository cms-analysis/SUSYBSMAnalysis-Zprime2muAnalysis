import FWCore.ParameterSet.Config as cms

process = cms.Process('RAWRECOSKIM')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:new.root'))
process.source.secondaryFileNames = cms.untracked.vstring('/store/data/Run2011A/SingleMu/RAW/v1/000/167/898/E468C1C7-7DA1-E011-97CB-BCAEC53296F4.root', '/store/data/Run2011A/SingleMu/RAW/v1/000/167/807/44F4E1B9-76A0-E011-8C58-0030487CD7E0.root')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.out = cms.OutputModule('PoolOutputModule',
                               fileName = cms.untracked.string('rawreco.root'),
                               outputCommands = cms.untracked.vstring('keep *', 'drop *_*_*_RAWRECOSKIM'),
                               )
process.outp = cms.EndPath(process.out)
