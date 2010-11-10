#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-120_7TeV-pythia6/dyzpforeff_dy120/9caa3d7638ff33984d7b458a78e3e8dd/pat_9_1_KwL.root']
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.options.wantSummary = True

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
tracks = ['global', 'inner', 'tpfms', 'picky', 'pmc', 'tmr', 'sigmaswitch']
rec_levels(process, tracks)

dy_gen_mass_cut = 'status == 3 && (pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039) && mass > %(lo)i && mass < %(hi)i'
process.DYGenMassFilter = cms.EDFilter('CandViewSelector',
                                       src = cms.InputTag('genParticles'),
                                       cut = cms.string(dy_gen_mass_cut % {'lo': 120, 'hi': 200}),
                                       filter = cms.bool(True),
                                       )

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')
process.EfficiencyFromMC.hardInteraction.src = cms.InputTag('genParticles')

process.p2 = cms.Path(process.Zprime2muAnalysisSequencePlain * process.EfficiencyFromMC)
process.p = cms.Path(process.DYGenMassFilter * process.Zprime2muAnalysisSequence)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi')

process.HistosFromPAT.leptonsFromDileptons = True
process.ResolutionUsingMC.leptonsFromDileptons = True
process.ResolutionUsingMC.hardInteraction.src = cms.InputTag('genParticles')
process.ResolutionUsingMC.doQoverP = True

process.p *= rec_level_module(process, process.HistosFromPAT,     'Histos', tracks)
process.p *= rec_level_module(process, process.ResolutionUsingMC, 'Resolution', tracks)

###############################################################################
    
import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(dataset)s
dbs_url = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
pset = histos_crab.py
total_number_of_events = -1
events_per_job = 50000

[USER]
ui_working_dir = crab/crab_ana_effres_%(name)s
return_data = 1
'''
    samples = [
        ('zmumu',  '/Zmumu/tucker-dyzpforeff_zmumu-9caa3d7638ff33984d7b458a78e3e8dd/USER', 20, 120),
        ('dy120',  '/DYToMuMu_M-120_7TeV-pythia6/tucker-dyzpforeff_dy120-9caa3d7638ff33984d7b458a78e3e8dd/USER',          120, 200),
        ('dy200',  '/DYToMuMu_M-200_7TeV-pythia6/tucker-dyzpforeff_dy200-9caa3d7638ff33984d7b458a78e3e8dd/USER',          200, 500),
        ('dy500',  '/DYToMuMu_M-500_7TeV-pythia6/tucker-dyzpforeff_dy500-9caa3d7638ff33984d7b458a78e3e8dd/USER',          500, 800),
        ('dy800',  '/DYToMuMu_M-800_7TeV-pythia6/tucker-dyzpforeff_dy800-9caa3d7638ff33984d7b458a78e3e8dd/USER',          800, 20000),
        ('zp1000', '/ZprimeSSMToMuMu_M-1000_7TeV-pythia6/tucker-dyzpforeff_zp1000-9caa3d7638ff33984d7b458a78e3e8dd/USER', -20000, 20000),
        ('zp1250', '/ZprimeSSMToMuMu_M-1250_7TeV-pythia6/tucker-dyzpforeff_zp1250-9caa3d7638ff33984d7b458a78e3e8dd/USER', -20000, 20000),
        ('zp1500', '/ZprimeSSMToMuMu_M-1500_7TeV-pythia6/tucker-dyzpforeff_zp1500-9caa3d7638ff33984d7b458a78e3e8dd/USER', -20000, 20000),
        ('zp1750', '/ZprimeSSMToMuMu_M-1750_7TeV-pythia6/tucker-dyzpforeff_zp1750-9caa3d7638ff33984d7b458a78e3e8dd/USER', -20000, 20000),
        ]

    just_testing = 'testing' in sys.argv

    if 'dump_files' in sys.argv:
        for sample in samples:
            print "('%s'," % sample[0]
            os.system('dbss ana02 find file where dataset=%s' % sample[1])
        sys.exit(0)

    for name, dataset, lo, hi in samples:
        open('crab.cfg', 'wt').write(crab_cfg % locals())

        new_py = open('histos.py').read()
        new_cut = dy_gen_mass_cut % locals()
        new_py += '\nprocess.DYGenMassFilter.cut = %(new_cut)s\n' % locals()
        open('histos_crab.py', 'wt').write(new_py)

        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm -v crab.cfg histos_crab.py')
