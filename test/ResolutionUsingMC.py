#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi')
process.p = cms.Path(process.Zprime2muAnalysisSequence * process.ResolutionUsingMC)
