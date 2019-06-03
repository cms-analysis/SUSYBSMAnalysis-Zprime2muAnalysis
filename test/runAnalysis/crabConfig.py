
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_pdf_CITo2Mu_Lam40TeVDesRR_M800to1300'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/MC17_CITo2Mu_M800to1300_CP5_Lam40TeVDesRRPythia8_v1/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_pdf_CITo2Mu_Lam40TeVDesRR_M800to1300'
config.Data.outLFNDirBase = '/store/user/jschulte/'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 8000
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

