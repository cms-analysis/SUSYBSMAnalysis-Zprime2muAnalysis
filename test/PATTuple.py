#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process

#from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import configure
#configure(process, muonsOnly=True) # useMonteCarlo=False)

process.source.fileNames = ['/store/relval/CMSSW_3_6_3/RelValZMM/GEN-SIM-RECO/START36_V10-v1/0006/369DA6CE-5878-DF11-8688-0030486790BE.root']
process.GlobalTag.globaltag = 'START36_V10::All'

process.p = cms.Path(process.patDefaultSequence)
