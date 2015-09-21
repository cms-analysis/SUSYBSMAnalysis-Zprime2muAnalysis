from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.workArea = ''
config.General.requestName = 'datamc_SingleMuonRun2015B-Prompt_251500_251642_20150717155341'
config.section_('JobType')
config.JobType.psetName = 'crab/psets/tuple_data_crab_SingleMuonRun2015B-Prompt_251500_251642_20150717155341.py'
config.JobType.pluginName = 'Analysis'
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/AOD'
config.Data.publishDBS = 'phys03'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/rradogna'
config.Data.splitting = 'LumiBased'
config.Data.inputDBS = 'global'
config.Data.lumiMask ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
config.Data.publishDataName = 'datamc_SingleMuonRun2015B-Prompt_251500_251642_20150717155341'
config.Data.publication = True
config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
config.section_('User')
config.section_('Debug')
