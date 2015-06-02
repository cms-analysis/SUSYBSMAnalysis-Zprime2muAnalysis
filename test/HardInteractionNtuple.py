#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
process.source.fileNames = ['/store/mc/Summer11/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/24D6F0D2-E77C-E011-8C27-003048D43944.root']

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
process.HardInteractionNtuple = cms.EDAnalyzer('HardInteractionNtuple', hardInteraction = hardInteraction, muon_src = cms.InputTag('muons'), picky_src = cms.InputTag('tevMuons', 'picky'))
process.HardInteractionNtuple.hardInteraction.src = 'genParticles'
process.p = cms.Path(process.HardInteractionNtuple)

if __name__ == '__main__' and 'submit' in sys.argv:
crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '%(name)s' 
config.General.workArea = 'HIntuple_%(name)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HardInteractionNtuple.py'   
config.JobType.priority = 1

config.Data.inputDataset = '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 10000
config.Data.publication = True
config.Data.publishDataName = '%(name)s'
config.Data.outLFN = '/store/user/federica/HIntuple' 

config.Site.storageSite = 'T2_US_Purdue'

'''
#### ??
    scheduler = 'condor'
    samples = [
        ('zmumu', '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy120', '/DYToMuMu_M-120_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy200', '/DYToMuMu_M-200_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ]

    just_testing = 'testing' in sys.argv

    for name, dataset in samples:
        open('crab.py', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab submit -c all')
            os.system('rm crab.py')
