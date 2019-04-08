
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_resolution_dy800to1400_2017'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_resolution_dy800to1400_2017'
config.Data.outLFNDirBase = '/store/user/zhangfa/ADD2016MC'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

