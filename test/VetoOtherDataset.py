#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.VetoOtherDataset = cms.EDFilter('VetoOtherDataset',
                                        hlt_process_name = cms.string('HLT'),
                                        dataset_to_veto = cms.string('MuMonitor'),
                                        )

process.pr = cms.EDAnalyzer('PrintEvent', trigger_results_src = cms.InputTag('TriggerResults', '', 'HLT'))
process.MessageLogger.categories.append('PrintEvent')

process.p = cms.Path(process.pr * process.VetoOtherDataset)

process.VetoOtherDataset2 = process.VetoOtherDataset.clone(dataset_to_veto='Mu')
process.p2 = cms.Path(process.VetoOtherDataset2)

process.VetoOtherDatasetMB = process.VetoOtherDataset.clone(dataset_to_veto='MinimumBias')
process.pMB = cms.Path(process.VetoOtherDatasetMB)

process.source.fileNames = ['file:crab/crab_datamc_Run2010A_DileptonMu/res/merged.root', 'file:crab/crab_datamc_Run2010B_DileptonMu/res/merged.root']
process.GlobalTag.globaltag = 'GR_R_38X_V15::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Nov4Run2010ABMuonsOnly
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*Nov4Run2010ABMuonsOnly)

