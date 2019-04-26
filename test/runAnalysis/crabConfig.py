
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_resolution_dyPt50To150_2Jets_2017_2017_whystuck'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/DY2JetsToLL_M-50_LHEZpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_resolution_dyPt50To150_2Jets_2017_whystuck'
config.Data.outLFNDirBase = '/store/user/zhangfa/DYMC2017resNOMUOrerun'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T3_US_FNALLPC'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

