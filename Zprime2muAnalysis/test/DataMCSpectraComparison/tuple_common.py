import os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process

process.p = cms.Path(process.countPatLeptons)

# Loose cut on muons; stronger cuts to be applied for different
# sets of plots (e.g. add our isolation cut, or apply VBTF).
process.selectedPatMuons.cut = 'isGlobalMuon && tunePMuonBestTrack().pt > 20'

process.countPatMuons.minNumber = 0
#process.countPatLeptons.electronSource = cms.InputTag('heepPatElectrons')
process.countPatLeptons.electronSource = cms.InputTag('patElectrons')
process.countPatLeptons.minNumber = 2

crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'datamc_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(pset)s'   
###config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
%(job_control)s
config.Data.publication = True
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/jschulte'

#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_DE_RWTH'
'''

os.system('mkdir -p crab/psets')
