import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.source.fileNames = ['/store/data/Run2015C/SingleMuon/AOD/PromptReco-v1/000/253/954/00000/E6980384-4441-E511-BB6E-02163E014728.root']
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v0'
process.maxEvents.input = 5000
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
process.CheckPrescale.dump_prescales = True

process.Mu27       = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu27_v2'))
#process.Mu15eta2p1 = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu15_eta2p1_v3', 'HLT_Mu15_eta2p1_v4'))
process.Mu24eta2p1 = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu24_eta2p1_v2'))

process.MessageLogger.suppressWarning = cms.untracked.vstring('Mu27', 'Mu24eta2p1')

#process.p = cms.Path(process.Mu17 * process.Mu15eta2p1 * process.Mu24eta2p1)
process.p = cms.Path(process.Mu24eta2p1*process.Mu27)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'getprescales_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'GetPrescales.py'
#config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 300
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_MuonPhys_v4.txt'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys.txt'
config.Data.publication = False
config.Data.publishDataName = 'getprescales_%(name)s'
config.Data.outLFNDirBase = '/store/user/rradogna'

#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_IT_Legnaro'
'''

    just_testing = 'testing' in sys.argv

    dataset_details = [
#        ('SingleMu2012A-13Jul2012', '/SingleMu/Run2012A-13Jul2012-v1/AOD'),
#        ('SingleMu2012B-13Jul2012', '/SingleMu/Run2012B-13Jul2012-v1/AOD'),
#        ('SingleMu2012C_24Aug2012', '/SingleMu/Run2012C-24Aug2012-v1/AOD'),
#        ('SingleMu2012C_Prompt',    '/SingleMu/Run2012C-PromptReco-v2/AOD'),
#        ('SingleMuon2015B_Prompt',    '/SingleMuon/Run2015B-PromptReco-v1/AOD'),
        ('SingleMuon2015C_Prompt',    '/SingleMuon/Run2015C-PromptReco-v1/AOD'),
        ]

    for name, dataset in dataset_details:
        print name
        open('crabConfig.py', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab submit -c crabConfig.py')

    if not just_testing:
        os.system('rm -f crabConfig.py tmp.json')
