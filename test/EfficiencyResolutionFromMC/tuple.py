#!/usr/bin/env python

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import cms, process
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName, AODOnly, pruneMCLeptons

process.maxEvents.input = 50
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
process.source.fileNames = ['/store/mc/RunIIFall15DR76/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/0094C520-D9B0-E511-B9FC-002590E3A024.root',
 #                           '/store/relval/CMSSW_7_4_0/RelValZpMM_2250_13TeV_Tauola/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/B6A4F83A-14DB-E411-8A01-0025905B8596.root'
]

process.p = cms.Path(process.countPatMuons)

pruneMCLeptons(process, use_sim=True)
AODOnly(process)

process.countPatMuons.minNumber = 0

process.out.outputCommands = [
    'drop *',
    'keep *_prunedMCLeptons_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep *_addPileupInfo_*_*',
    'keep patMuons_cleanPatMuonsTriggerMatch__PAT',
    'keep L1GlobalTriggerObjectMaps_l1L1GtObjectMap_*_*',
    'keep L1GlobalTriggerReadoutRecord_gtDigis__RECO',
    'keep triggerTriggerEvent_hltTriggerSummaryAOD__HLT*',
    'keep triggerTriggerEvent_hltTriggerSummaryAOD__REDIGI*',
    'keep patElectrons_cleanPatElectrons__PAT',
    'keep patPhotons_cleanPatPhotons__PAT',
    'keep edmTriggerResults_TriggerResults__HLT*',
    'keep edmTriggerResults_TriggerResults__REDIGI*',
    'keep GenEventInfoProduct_generator__SIM',
    #'keep GenEventInfoProduct_generator__HLT',
    'keep edmTriggerResults_TriggerResults__PAT',
    'keep *_patTrigger_*_*', # keep these two for now, for Level-1 decisions
    'keep *_patTriggerEvent_*_*',
    ]

import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'effres_%(name)s' 
config.General.workArea = 'PAT_effres_%(name)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(pset_fn)s'   
###config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 10000
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'effres_%(name)s'
config.Data.outLFNDirBase = '/store/user/alfloren/PAATuples'
config.Data.ignoreLocality = True 

config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ["T2_CH_CERN"]

'''
    #os.system('mkdir -p psets crab')
    
    just_testing = 'testing' in sys.argv
    
    samples = [
       # ('dy50to120',   '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
       # ('dy120to200',  '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
       # ('dy200to400',  '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
        ('dy400to800',  '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
        #('dy800to1400',  '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
       # ('dy1400to2300', '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
       # ('dy2300to3500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
       # ('dy3500to4500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'),
        #('dy4500to6000', '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/AODSIM'),

  ]

    for name, dataset in samples:
        print name
        events_per_job = 5000

        pset = open('tuple.py').read()
        #if name == 'dy20':
        #    pset += '\nswitchHLTProcessName(process, "REDIGI38X")\n'
        pset_fn = 'tuple_effres_crab_%s.py' % name
        open(pset_fn, 'wt').write(pset)
        
        open('crabConfig.py', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab submit -c crabConfig.py')
            # os.system('rm crab.py')
