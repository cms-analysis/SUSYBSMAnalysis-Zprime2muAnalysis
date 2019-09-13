
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons_dy6000toInf'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ZToEE_NNPDF31_13TeV-powheg_M_6000_Inf/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_electrons_dy6000toInf'
config.Data.outLFNDirBase = '/store/user/jschulte/'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
#config.General.instance = 'preprod' 
config.Site.whitelist = ["T2_US_*"]
config.Site.blacklist = ['T2_US_Caltech']
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 8000
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

