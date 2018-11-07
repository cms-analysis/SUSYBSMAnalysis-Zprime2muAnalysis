
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_resolution_dy_3Jets_2017'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_resolution_dy_3Jets_2017'
config.Data.outLFNDirBase = '/store/user/jschulte'
#config.Data.ignoreLocality = True 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 400000

