#!/usr/bin/env python

files = [
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/EAB4AE8B-FEC8-DF11-B158-0018FE284DE8.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/E870712C-D2C8-DF11-867E-00163E061001.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/CCD6E542-DFC8-DF11-933E-00215E2EB700.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/A07DC8B0-D2C8-DF11-95C4-00237DF29408.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/A000027B-D3C8-DF11-A571-00163E021401.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/9A7B86B0-2BC9-DF11-B67B-0018FE286F12.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/98E52639-27C9-DF11-87A0-0018FE284D7C.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/901DF82D-D2C8-DF11-890D-00163E071301.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/8C2A9023-95C9-DF11-8769-0018FE283D46.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/66CF2604-9AC9-DF11-88BB-00237DF2B4E0.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/6469CFAD-00C9-DF11-B8CD-0018FE284C14.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/14B530EB-92C9-DF11-8A27-00163E031401.root',
    '/store/mc/Fall10/DYToMuMu_M-200_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0002/12EF6BBC-2BC9-DF11-898D-0018FE283CCC.root',
    ]
maxEvents = -1
output_filename = 'Parametrizer.root'

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
