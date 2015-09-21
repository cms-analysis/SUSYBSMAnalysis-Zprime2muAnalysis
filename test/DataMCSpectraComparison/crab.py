

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_datamc_Run2012MuonsOnly_SingleMuRun2012A_13Jul2012_190450_193751'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
#config.JobType.priority = 1

config.Data.inputDataset =  '/SingleMu/slava-datamc_SingleMuRun2012A-13Jul2012_190450_193751_20121011073628-426a2d966f78bce6bde85f3ed41c07ba/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased' 

total_number_of_lumis = -1
lumis_per_job = 500
lumi_mask = tmp.json
config.Data.publication = False
config.Data.publishDataName = 'ana_datamc_Run2012MuonsOnly_SingleMuRun2012A_13Jul2012_190450_193751'
config.Data.outLFNDirBase = '/store/user/rradogna'

#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_IT_Legnaro'

