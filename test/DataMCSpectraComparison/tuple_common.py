import os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.p = cms.Path(process.countPatLeptons)

# Loose cut on muons; stronger cuts to be applied for different
# sets of plots (e.g. add our isolation cut, or apply VBTF).
process.selectedPatMuons.cut = 'isGlobalMuon && pt > 20'

# Want to select only events that have at least two leptons (=
# muons+electrons), where the electrons must pass HEEP id, but don't
# want to force HEEP id on selectedPatElectrons so as not to screw up
# the jet cleaning until we study this.
process.heepPatElectrons = cms.EDFilter('PATElectronSelector',
                                        src = cms.InputTag('patElectrons'),
                                        cut = cms.string('userInt("HEEPId") == 0')
                                        )

process.patDefaultSequence.replace(process.selectedPatElectrons, process.selectedPatElectrons * process.heepPatElectrons)
process.countPatMuons.minNumber = 0
process.countPatLeptons.electronSource = cms.InputTag('heepPatElectrons')
process.countPatLeptons.minNumber = 2

crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '%(name)s' 
config.General.workArea = 'PAT_%(name)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(pset)s'   
config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 10000
config.Data.publication = True
config.Data.publishDataName = '%(name)s'
config.Data.outLFNDirBase = '/store/user/federica/PATTuple' 

config.Site.storageSite = 'T2_US_Purdue'

'''

#os.system('mkdir -p crab/psets')
