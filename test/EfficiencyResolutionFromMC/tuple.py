#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.maxEvents.input = 100
process.GlobalTag.globaltag = 'START37_V6::All'
process.source.fileNames = ['/store/mc/Spring10/DYToMuMu_M-120_7TeV-pythia6/GEN-SIM-RECO/START3X_V26-v2/0026/F44E2CAC-A35A-DF11-B61B-002481CFE804.root']
process.p = cms.Path(process.patDefaultSequence)

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run36xOn35xInput
run36xOn35xInput(process, 'ak5GenJets')

process.countPatMuons.minNumber = 0

process.out.outputCommands += [
    'drop *_patTrigger_*_*',
    'drop *_patTriggerEvent_*_*',
    'drop *_prunedGenSimLeptons_*_*',
    'keep *_genParticles_*_*',
    'drop *_cleanPatTaus_*_*'
    ]

import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = %(scheduler)s

[CMSSW]
datasetpath = %(dataset)s
pset = tuple.py
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

[GRID]
#ce_black_list = ufl
'''

    just_testing = 'testing' in sys.argv
    
    samples = [
        ('zmumu', '/Zmumu/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO'),
        ('dy120', '/DYToMuMu_M-120_7TeV-pythia6/Spring10-START3X_V26-v2/GEN-SIM-RECO'),
        ('dy200', '/DYToMuMu_M-200_7TeV-pythia6/Spring10-START3X_V26-v2/GEN-SIM-RECO'),
        ('dy500', '/DYToMuMu_M-500_7TeV-pythia6/Spring10-START3X_V26-v2/GEN-SIM-RECO'),
        ('dy800', '/DYToMuMu_M-800_7TeV-pythia6/Spring10-START3X_V26-v2/GEN-SIM-RECO'),
        ('zp1000', '/ZprimeSSMToMuMu_M-1000_7TeV-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO'),
        ('zp1250', '/ZprimeSSMToMuMu_M-1250_7TeV-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO'),
        ('zp1500', '/ZprimeSSMToMuMu_M-1500_7TeV-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO'),
        ('zp1750', '/ZprimeSSMToMuMu_M-1750_7TeV-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO'),
        ]

    scheduler = 'condor'
    for name, dataset in samples:
        events = 5000 if name != 'zmumu' else 20000
        #scheduler = 'condor' if name == 'zmumu' else 'glite'

        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -cfg %s -create -submit all' % crab_fn)
    if not just_testing:
        os.system('rm -v crab.cfg')
