
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_muons_2016_SingleMuonRun2016B-23Sep2016_v3'
config.General.workArea = 'crabRecovery'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = 'dileptonAna_muons_2016_SingleMuonRun2016B-23Sep2016_v3'
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
config.Data.lumiMask = '/afs/cern.ch/work/j/jschulte/test/CMSSW_10_2_15_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/runAnalysis/crab/crab_dileptonAna_muons_2016_SingleMuonRun2016B-23Sep2016_v3/results/notFinishedLumis.json'

