
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_genMassOther_electrons_tW'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = 'dileptonAna_genMassOther_electrons_tW'
config.Data.outLFNDirBase = '/store/user/jschulte/'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T2_US_Purdue'
#config.JobType.maxMemoryMB  = 4000
config.JobType.allowUndistributedCMSSW = True
config.Site.blacklist = ['T2_US_Caltech']

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

