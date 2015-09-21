#!/usr/bin/env python

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import cms, process
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName, AODOnly, pruneMCLeptons

process.maxEvents.input = 50
process.GlobalTag.globaltag = 'MCRUN2_74_V9A'
process.source.fileNames = ['/store/relval/CMSSW_7_4_0/RelValZpMM_2250_13TeV_Tauola/GEN-SIM-RECO/MCRUN2_74_V7-v1/00000/AE8D58C2-14DB-E411-A038-002618943901.root',
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
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(pset_fn)s'
###config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = %(events_per_job)s
config.Data.publication = True
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'effres_%(name)s'
config.Data.outLFNDirBase = '/store/user/rradogna'

#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_IT_Legnaro'
'''
    #os.system('mkdir -p psets crab')
    
    just_testing = 'testing' in sys.argv
    
    samples = [
        ('dy50',   '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM'),
#        ('dy120',  '/DYToMuMu_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy200',  '/DYToMuMu_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy500',  '/DYToMuMu_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy800',  '/DYToMuMu_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy1000', '/DYToMuMu_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy1500', '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy2000', '/DYToMuMu_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('dy9500',  '/DYJetsToEEMuMu_M-9500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM'),
#        ('dy200_c1',  '/DYToMuMu_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('dy500_c1',  '/DYToMuMu_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('dy800_c1',  '/DYToMuMu_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('dy1000_c1', '/DYToMuMu_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('dy1500_c1', '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('dy2000_c1', '/DYToMuMu_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('dy120_c2',  '/DYToMuMu_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('dy200_c2',  '/DYToMuMu_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('dy500_c2',  '/DYToMuMu_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('dy800_c2',  '/DYToMuMu_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('dy1000_c2', '/DYToMuMu_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('dy1500_c2', '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('dy2000_c2', '/DYToMuMu_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'),
#        ('zp750',  '/ZprimePSIToMuMu_M-750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp1000', '/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp1250', '/ZprimePSIToMuMu_M-1250_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp1500', '/ZprimePSIToMuMu_M-1500_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp1750', '/ZprimePSIToMuMu_M-1750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp2000', '/ZprimePSIToMuMu_M-2000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp2250', '/ZprimePSIToMuMu_M-2250_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp2500', '/ZprimePSIToMuMu_M-2500_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp2750', '/ZprimePSIToMuMu_M-2750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp3000', '/ZprimePSIToMuMu_M-3000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'),
#        ('zp5000',  '/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM'),
#        ('zp1000_c1', '/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp1250_c1', '/ZprimePSIToMuMu_M-1250_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp1500_c1', '/ZprimePSIToMuMu_M-1500_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp1750_c1', '/ZprimePSIToMuMu_M-1750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp2000_c1', '/ZprimePSIToMuMu_M-2000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp2250_c1', '/ZprimePSIToMuMu_M-2250_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp2500_c1', '/ZprimePSIToMuMu_M-2500_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp2750_c1', '/ZprimePSIToMuMu_M-2750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
#        ('zp3000_c1', '/ZprimePSIToMuMu_M-3000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM'),
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
            #os.system('rm crabConfig.py')
