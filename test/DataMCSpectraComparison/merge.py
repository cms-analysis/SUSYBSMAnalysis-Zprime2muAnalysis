#!/usr/bin/env python

import os, re, sys, glob, FWCore.ParameterSet.Config as cms
from pprint import pprint

process = cms.Process('MergePAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Look for just a list of files in argv first.
files = ['file:%s' % x for x in sys.argv if os.path.isfile(x) and '.root' in x]
name = 'merged.root'
if not files:
    # Else, files from crab dir mode.
    from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import last_crab_dir, files_from_crab_dir
    crab_dir = [x for x in sys.argv if os.path.isdir(x)][0]
    if not crab_dir:
        crab_dir = last_crab_dir()
        print 'Auto-crab dir:', crab_dir
    files = files_from_crab_dir(crab_dir)
    name = os.path.join(crab_dir, 'res', 'merged.root')
    
print 'Files to run over:'
pprint(files)
process.source.fileNames = files
print 'Merging to', name
process.out = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string(name))
process.outp = cms.EndPath(process.out)
