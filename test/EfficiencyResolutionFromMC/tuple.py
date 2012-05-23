#!/usr/bin/env python

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import cms, process
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName, AODOnly

process.maxEvents.input = 50
process.GlobalTag.globaltag = 'START52_V9::All'
process.source.fileNames = ['file:/uscms/home/tucker/nobackup/store/mc/Summer11/ZprimeSSMToMuMu_M-2250_TuneZ2_7TeV-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/7032435D-1A93-E011-94D5-0017A4770008.root']
process.p = cms.Path(process.patDefaultSequence)

AODOnly(process)

process.countPatMuons.minNumber = 0

process.out.outputCommands = [
    'drop *',
    'keep *_prunedMCLeptons_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep *_addPileupInfo_*_*',
    'keep patMuons_cleanPatMuonsTriggerMatch__PAT',
    'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT*',
    'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__REDIGI*',
    'keep L1GlobalTriggerReadoutRecord_gtDigis__RECO',
    'keep triggerTriggerEvent_hltTriggerSummaryAOD__HLT*',
    'keep triggerTriggerEvent_hltTriggerSummaryAOD__REDIGI*',
    'keep patElectrons_cleanPatElectrons__PAT',
    'keep patPhotons_cleanPatPhotons__PAT',
    'keep edmTriggerResults_TriggerResults__HLT*',
    'keep edmTriggerResults_TriggerResults__REDIGI*',
    'keep GenEventInfoProduct_generator__HLT',
    'keep edmTriggerResults_TriggerResults__PAT',
    ]

import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(dataset)s
pset = %(pset_fn)s
total_number_of_events = -1
events_per_job = %(events_per_job)s
get_edm_output = 1

[USER]
ui_working_dir = crab/crab_effres_%(name)s
copy_data = 1
storage_element = T3_US_FNALLPC
check_user_remote_dir = 0
publish_data = 1
publish_data_name = effres_%(name)s
dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
'''

    os.system('mkdir -p psets crab')
    
    just_testing = 'testing' in sys.argv
    
    samples = [
        ('dy20',   '/DYToMuMu_M_20_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy120',  '/DYToMuMu_M_120_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy200',  '/DYToMuMu_M_200_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy500',  '/DYToMuMu_M_500_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy800',  '/DYToMuMu_M_800_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy1000', '/DYToMuMu_M_1000_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy1300', '/DYToMuMu_M-1300_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('dy1600', '/DYToMuMu_M-1600_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp750',  '/ZprimePSIToMuMu_M-750_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp1000', '/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp1250', '/ZprimePSIToMuMu_M-1250_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp1500', '/ZprimePSIToMuMu_M-1500_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp1750', '/ZprimePSIToMuMu_M-1750_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp2000', '/ZprimePSIToMuMu_M-2000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp2250', '/ZprimePSIToMuMu_M-2250_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp2500', '/ZprimePSIToMuMu_M-2500_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp2750', '/ZprimePSIToMuMu_M-2750_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        ('zp3000', '/ZprimePSIToMuMu_M-3000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM'),
        # RSG samples...
    ]

    for name, dataset in samples:
        print name
        events_per_job = 26000

        pset = open('tuple.py').read()
        #if name == 'dy20':
        #    pset += '\nswitchHLTProcessName(process, "REDIGI38X")\n'
        pset_fn = 'psets/tuple_effres_crab_%s.py' % name
        open(pset_fn, 'wt').write(pset)
        
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg')
