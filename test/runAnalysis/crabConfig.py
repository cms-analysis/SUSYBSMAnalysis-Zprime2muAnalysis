
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons_DoubleEG2017F-31Mar2018-v1'
config.General.workArea = 'crabNew3'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = 'dileptonAna_electrons_DoubleEG2017F-31Mar2018-v1'
config.Data.outLFNDirBase = '/store/user/jschulte/'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 4000
config.JobType.allowUndistributedCMSSW = True
config.Site.blacklist = ['T2_US_Caltech']

config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'

