#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V13::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import files_from_argv
files_from_argv(process)

process.load('HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi')
process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')

process.p = cms.Path(process.hltEventAnalyzerAOD * process.triggerSummaryAnalyzerAOD)

#from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
#switch_hlt_process_name(process, 'duh')
