
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons__2016DoublEG2016H-17Jul2017-v1'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/DoubleEG/Run2016H-17Jul2018-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_electrons__2016DoublEG2016H-17Jul2017-v1'
config.Data.outLFNDirBase = '/store/user/jschulte'
#config.Data.ignoreLocality = True 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_US_Purdue'
#config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 200
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

