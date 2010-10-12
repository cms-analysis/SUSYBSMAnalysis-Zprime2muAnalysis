#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler',
                                       hlt_src = cms.InputTag('TriggerResults', '', 'HLT'),
                                       dimu_src = cms.InputTag('dimuons')
                                       )

process.p = cms.Path(process.Zprime2muAnalysisSequence * process.SimpleNtupler)

from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import files_from_argv
files_from_argv(process)

lumis_to_process = open('/afs/cern.ch/user/t/tucker/runreg/output/full.cmssw').read().split(',')
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*lumis_to_process)
