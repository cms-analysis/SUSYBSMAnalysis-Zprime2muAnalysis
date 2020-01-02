
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons_2018_CITo2E_Lam100kTeVDesRL_M1300to2000'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/MC17_CITo2E_M1300to2000_CP5_Lam100kTeVDesLRPythia8_v1/RunIIFall17MiniAOD-PU2017_94X_mc2017_realistic_v11-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_electrons_2018_CITo2E_Lam100kTeVDesRL_M1300to2000'
config.Data.outLFNDirBase = '/store/user/jschulte/'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
#config.General.instance = 'preprod' 
config.Site.whitelist = ["T2_US_*"]
config.Site.blacklist = ['T2_US_Caltech']
config.Site.storageSite = 'T2_US_Purdue'
#config.JobType.maxMemoryMB  = 4000
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

