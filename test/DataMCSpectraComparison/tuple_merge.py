#!/usr/bin/env python

import datetime, sys, os, FWCore.ParameterSet.Config as cms
from collections import defaultdict

process = cms.Process('MergePAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:pat.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 30000

process.out = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('merged.root'))
process.outp = cms.EndPath(process.out)

process.options.emptyRunLumiMode = cms.untracked.string('doNotHandleEmptyRunsAndLumis')
process.source.inputCommands = cms.untracked.vstring('keep *', 'drop *_MEtoEDMConverter_*_*')
process.out.outputCommands = cms.untracked.vstring('keep *', 'drop LumiDetails_lumiProducer_*_*', 'drop LumiSummary_lumiProducer_*_*', 'drop RunSummary_lumiProducer_*_*')

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'merge_%(uniq)s_%(name)s' 
config.General.workArea = 'PAT_merge_%(uniq)s_%(name)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'tuple_merge.py'   
config.JobType.priority = 1

config.Data.inputDataset = '%(ana_dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 10000
config.Data.publication = True
config.Data.publishDataName = 'merge_%(uniq)s_%(name)s'
config.Data.outLFN = '/store/user/federica/PATTuple' 

config.Site.storageSite = 'T2_US_Purdue'

'''

    just_testing = 'testing' in sys.argv
    uniq = datetime.datetime.today().strftime('%Y%m%d%H%M%S')

    cfgs = []
    if True:
        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import wjets, inclmu15
        to_merge = [wjets, inclmu15]
        for sample in to_merge:
            sample.uniq = uniq
            sample.job_control = 'total_number_of_events = -1'
            cfgs.append(crab_cfg % sample)
    else:
        to_merge = [
#            ('prompt165071165558', '/SingleMu/tucker-datamc_SingleMu2011A_prompt_165071_165558_20110526195544-8788f1b70631d1fb57e97a89f5e8007c/USER'),
#            ('prompt165559165627', '/SingleMu/tucker-datamc_SingleMu2011A_prompt_165559_165627_20110528141442-8788f1b70631d1fb57e97a89f5e8007c/USER'),
            ('may10temp',          '/SingleMu/tucker-datamc_SingleMuRun2011A_May10-7b18eaf160f3796dde5c7353a5ebfddf/USER'),
            ]
        for name, ana_dataset in to_merge:
            job_control = 'total_number_of_lumis = -1'
            if name == 'may10temp':
                job_control += '\nlumi_mask = may10tempdrop.json'
            cfgs.append(crab_cfg % locals())

    for cfg in cfgs:
        open('crab.py', 'wt').write(cfg)
        if not just_testing:
            os.system('crab submit -c all')
        else:
            print cfg

    if not just_testing:
        os.system('rm crab.py')
