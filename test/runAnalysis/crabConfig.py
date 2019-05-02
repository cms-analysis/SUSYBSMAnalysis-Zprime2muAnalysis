
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_muons_dy50to120_2017_whystuck'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_muons_dy50to120_whystuck'
config.Data.outLFNDirBase = '/store/user/zhangfa/DYMC2017resNOMUOrerun'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T3_US_FNALLPC'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

