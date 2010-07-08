#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import configure
configure(process, useMonteCarlo=False)

process.source.fileNames = ['file:csonia.root']
process.GlobalTag.globaltag = 'GR10_P_V7::All'

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
#process.selectedPatMuons.cut = 'muonID("GlobalMuonPromptTight") && muonID("TMOneStationLoose") && (globalTrack.hitPattern.numberOfValidMuonCSCHits + globalTrack.hitPattern.numberOfValidMuonDTHits) >= 1 && innerTrack.hitPattern.trackerLayersWithMeasurement >= 6'

process.p = cms.Path(process.goodData * process.patDefaultSequence)
