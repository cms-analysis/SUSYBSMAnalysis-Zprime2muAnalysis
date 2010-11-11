#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import files_from_argv
files_from_argv(process)

hlt = 'REDIGI38X'

process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD', '', hlt)

process.printEvent = cms.EDAnalyzer('PrintEvent', hlt_src = cms.InputTag('TriggerResults', '', hlt))
process.MessageLogger.categories.append('PrintEvent')

process.p = cms.Path(process.printEvent * process.triggerSummaryAnalyzerAOD)
