#!/usr/bin/env python

import os, re, sys, glob, FWCore.ParameterSet.Config as cms
from pprint import pprint

process = cms.Process('MergePAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import last_crab_dir, files_from_crab_dir

crab_dir = [x for x in sys.argv if os.path.isdir(x)][0]
if not crab_dir:
    crab_dir = last_crab_dir()

files = files_from_crab_dir(crab_dir)
print 'Files for %s to run over:' % crab_dir
pprint(files)
process.source.fileNames = files

name = os.path.join(crab_dir, 'res', 'merged.root')
print 'Merging to', name
process.out = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string(name))
process.outp = cms.EndPath(process.out)
