
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons_2018_DoubleEGRun2018D-PromptReco-v2'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/EGamma/Run2018D-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_electrons_2018_DoubleEGRun2018D-PromptReco-v2'
config.Data.outLFNDirBase = '/store/user/jschulte'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 400
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

