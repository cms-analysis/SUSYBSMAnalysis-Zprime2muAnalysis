import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.source.fileNames = ['/store/data/Run2012C/SingleMu/AOD/PromptReco-v2/000/202/272/94E42528-89F9-E111-BE44-BCAEC53296F8.root']
process.GlobalTag.globaltag = 'GR_P_V42_AN2::All'
process.maxEvents.input = 5000
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
process.CheckPrescale.dump_prescales = True

process.Mu17       = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu17_v3'))
process.Mu15eta2p1 = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu15_eta2p1_v3', 'HLT_Mu15_eta2p1_v4'))
process.Mu24eta2p1 = process.CheckPrescale.clone(trigger_paths=cms.vstring('HLT_Mu24_eta2p1_v3', 'HLT_Mu24_eta2p1_v4', 'HLT_Mu24_eta2p1_v5'))

process.MessageLogger.suppressWarning = cms.untracked.vstring('Mu17', 'Mu15eta2p1', 'Mu24eta2p1')

#process.p = cms.Path(process.Mu17 * process.Mu15eta2p1 * process.Mu24eta2p1)
process.p = cms.Path(process.Mu24eta2p1)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '%(name)s' 
config.General.workArea = 'Prescale_%(name)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'GetPrescales.py'   
config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 10000
config.Data.publication = True
config.Data.publishDataName = '%(name)s'
config.Data.outLFN = '/store/user/federica/Prescale' 

config.Site.storageSite = 'T2_US_Purdue'

'''

    just_testing = 'testing' in sys.argv

    dataset_details = [
        ('SingleMu2012A-13Jul2012', '/SingleMu/Run2012A-13Jul2012-v1/AOD'),
        ('SingleMu2012B-13Jul2012', '/SingleMu/Run2012B-13Jul2012-v1/AOD'),
        ('SingleMu2012C_24Aug2012', '/SingleMu/Run2012C-24Aug2012-v1/AOD'),
        ('SingleMu2012C_Prompt',    '/SingleMu/Run2012C-PromptReco-v2/AOD'),
        ('SingleMu2012D_Prompt',    '/SingleMu/Run2012D-PromptReco-v1/AOD'),
        ]

    for name, dataset in dataset_details:
        print name
        open('crab.py', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab submit -c all')

    if not just_testing:
        os.system('rm -f crab.py tmp.json')
