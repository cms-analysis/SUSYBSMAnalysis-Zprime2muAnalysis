#!/usr/bin/env python

files = ['file:/uscms_data/d1/tucker/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN.root']
maxEvents = -1
output_filename = 'AsymmetryParametrizer.root'

################################################################################

import sys, FWCore.ParameterSet.Config as cms
process = cms.Process('AsymmetryParametrizer')

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.AsymmetryParametrizer_cfi')
process.p = cms.Path(process.AsymmetryParametrizer)

if 'assemble' in sys.argv:
    process.source = cms.Source('EmptySource')
    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
    process.AsymmetryParametrizer.assemble_only = True
else:
    process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring(*files))
    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))
    process.TFileService = cms.Service('TFileService', fileName = cms.string(output_filename))

    process.load('FWCore.MessageLogger.MessageLogger_cfi')
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
