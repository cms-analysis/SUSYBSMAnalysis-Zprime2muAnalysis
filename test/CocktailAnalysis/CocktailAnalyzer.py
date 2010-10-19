#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import files_from_argv

process.MuonsFromDimuons = cms.EDProducer('MuonsFromDimuons', dimuon_src = cms.InputTag('dimuons'))
process.CocktailAnalyzerAll     = cms.EDAnalyzer('CocktailAnalyzer', muon_src = cms.InputTag('MuonsFromDimuons'), eta_min = cms.double(0.0), eta_max = cms.double(3.0))
process.CocktailAnalyzerBarrel  = cms.EDAnalyzer('CocktailAnalyzer', muon_src = cms.InputTag('MuonsFromDimuons'), eta_min = cms.double(0.0), eta_max = cms.double(0.8))
process.CocktailAnalyzerOverlap = cms.EDAnalyzer('CocktailAnalyzer', muon_src = cms.InputTag('MuonsFromDimuons'), eta_min = cms.double(0.8), eta_max = cms.double(1.2))
process.CocktailAnalyzerEndcap  = cms.EDAnalyzer('CocktailAnalyzer', muon_src = cms.InputTag('MuonsFromDimuons'), eta_min = cms.double(1.2), eta_max = cms.double(3.0))

process.p = cms.Path(process.goodDataFilter * process.Zprime2muAnalysisSequence * process.MuonsFromDimuons * process.CocktailAnalyzerAll * process.CocktailAnalyzerBarrel * process.CocktailAnalyzerOverlap * process.CocktailAnalyzerEndcap)

files_from_argv(process)
process.GlobalTag.globaltag = 'GR10_P_V7::All'
process.TFileService.fileName = 'cocktail.root'
process.options.wantSummary = True
