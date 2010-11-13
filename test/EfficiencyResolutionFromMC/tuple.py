#!/usr/bin/env python

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import cms, process
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName

process.maxEvents.input = 50
process.GlobalTag.globaltag = 'START38_V12::All'
process.source.fileNames = ['/store/mc/Fall10/ZprimeSSMToMuMu_M-750_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0001/F034C4E4-4BCD-DF11-B652-00E0812EEEC3.root']
process.p = cms.Path(process.patDefaultSequence)

process.countPatMuons.minNumber = 0

process.out.outputCommands = [
    'drop *',
    'keep *_prunedGenSimLeptons_*_*',
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
events_per_job = %(events)s
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
        ('dy20',   '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('dy120',  '/DYToMuMu_M-120_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('dy200',  '/DYToMuMu_M-200_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('dy500',  '/DYToMuMu_M-500_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('dy800',  '/DYToMuMu_M-800_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('zp500',  '/ZprimeSSMToMuMu_M-500_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('zp750',  '/ZprimeSSMToMuMu_M-750_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('zp1000', '/ZprimeSSMToMuMu_M-1000_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('zp1250', '/ZprimeSSMToMuMu_M-1250_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('zp1500', '/ZprimeSSMToMuMu_M-1500_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('zp1750', '/ZprimeSSMToMuMu_M-1750_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
    ]

    for name, dataset in samples:
        if name == 'dy20':
            continue
        
        events = 11000 if name != 'dy20' else 40000
        pset = open('tuple.py').read()
        if name == 'dy20':
            pset += '\nswitchHLTProcessName(process, "REDIGI38X")\n'
        pset_fn = 'psets/tuple_effres_crab_%s.py' % name
        open(pset_fn, 'wt').write(pset)
        
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg')
