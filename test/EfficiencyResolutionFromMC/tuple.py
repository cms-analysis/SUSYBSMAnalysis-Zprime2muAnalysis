#!/usr/bin/env python

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import cms, process
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName, AODOnly

process.maxEvents.input = 50
process.GlobalTag.globaltag = 'START42_V11::All'
process.source.fileNames = ['file:/uscms/home/tucker/nobackup/store/mc/Summer11/ZprimeSSMToMuMu_M-1000_TuneZ2_7TeV-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/065F848B-3D93-E011-8E9F-78E7D164BFC8.root']
process.p = cms.Path(process.patDefaultSequence)

AODOnly(process)

process.countPatMuons.minNumber = 0

process.out.outputCommands = [
    'drop *',
    'keep *_prunedMCLeptons_*_*',
    'keep patMuons_cleanPatMuonsTriggerMatch__PAT',
    'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT*',
    'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__REDIGI*',
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
        ('dy20',   '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy120',  '/DYToMuMu_M-120_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy200',  '/DYToMuMu_M-200_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy500',  '/DYToMuMu_M-500_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy800',  '/DYToMuMu_M-800_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('dy1000', '/DYToMuMu_M-1000_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM'),
        ('zp750',  '/ZprimeSSMToMuMu_M-750_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        ('zp1000', '/ZprimeSSMToMuMu_M-1000_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        ('zp1250', '/ZprimeSSMToMuMu_M-1250_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        ('zp1500', '/ZprimeSSMToMuMu_M-1500_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        ('zp1750', '/ZprimeSSMToMuMu_M-1750_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        ('zp2000', '/ZprimeSSMToMuMu_M-2000_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        ('zp2250', '/ZprimeSSMToMuMu_M-2250_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM'),
        # RSG samples...
    ]

    for name, dataset in samples:
        events_per_job = 22000

        pset = open('tuple.py').read()
        #if name == 'dy20':
        #    pset += '\nswitchHLTProcessName(process, "REDIGI38X")\n'
        pset_fn = 'psets/tuple_effres_crab_%s.py' % name
        open(pset_fn, 'wt').write(pset)
        
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg')
