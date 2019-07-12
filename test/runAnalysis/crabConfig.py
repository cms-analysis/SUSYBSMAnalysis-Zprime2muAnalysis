
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons_GJets_HT600ToInf'
config.General.workArea = 'crabNew'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_electrons_GJets_HT600ToInf'
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

