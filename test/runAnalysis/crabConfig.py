
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_muons_2016_ADDGravToLL_LambdaT100k_M1700'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ADDGravToLL_LambdaT-100kTeV_M-1700ToInf_13TeV-pythia8/jschulte-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1-56fe694eb6687338c7aeb0151d5792b3/USER'
config.Data.inputDBS = 'phys03'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_muons_2016_ADDGravToLL_LambdaT100k_M1700'
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

