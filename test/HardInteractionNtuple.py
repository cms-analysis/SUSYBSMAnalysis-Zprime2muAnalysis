#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
process.source.fileNames = ['file:/uscms/home/tucker/scratch/zmumu.root']

from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
process.HardInteractionNtuple = cms.EDAnalyzer('HardInteractionNtuple', hardInteraction = hardInteraction)
process.HardInteractionNtuple.hardInteraction.src = 'genParticles'
process.p = cms.Path(process.HardInteractionNtuple)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = %(scheduler)s

[CMSSW]
datasetpath = %(dataset)s
pset = HardInteractionNtuple.py
total_number_of_events = -1
events_per_job = 50000

[USER]
ui_working_dir = crab/crab_genntuple_%(name)s
return_data = 1
'''

    samples = [
        ('pythia36',     '/Zmumu/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO'),
        ('powheg36',     '/DYToMuMu_M-20_7TeV-powheg-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO'),
        ('pythia',       '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('powheg',       '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall10-START38_V12-v1/GEN-SIM-RECO'),
        ('powheg36ct66', '/Zmumu_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v2/GEN-SIM-RECO'),
        ]

    just_testing = 'testing' in sys.argv

    for name, dataset in samples:
        scheduler = 'glite' if name == 'pythia36' else 'condor'
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg')
