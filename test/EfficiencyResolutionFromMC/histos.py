#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process
process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/effres_dy20/8ca75260210b8943d361f4da5b0c0bcc/pat_1_1_Zi3.root']
process.options.wantSummary = True

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
tracks = ['global', 'inner', 'tpfms', 'picky', 'pmc', 'tmr', 'sigmaswitch']
rec_levels(process, tracks)

dy_gen_mass_cut = 'status == 3 && (pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039) && mass > %(lo)i && mass < %(hi)i'
process.DYGenMassFilter = cms.EDFilter('CandViewSelector',
                                       src = cms.InputTag('prunedGenSimLeptons'),
                                       cut = cms.string(dy_gen_mass_cut % {'lo': 0, 'hi': 10000}),
                                       filter = cms.bool(True),
                                       )

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')

process.HLTSingleObjects = cms.EDProducer('HLTLeptonsFromTriggerEvent',
                                          summary = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
                                          leptons = cms.VInputTag(cms.InputTag('hltL3MuonCandidates', '', 'HLT'))
                                          )
process.EfficiencyFromMC.hlt_obj_src = 'HLTSingleObjects'
process.EfficiencyFromMC.hlt_single_min_pt = 15

process.p2 = cms.Path(process.Zprime2muAnalysisSequencePlain * process.HLTSingleObjects * process.EfficiencyFromMC)
process.p = cms.Path(process.DYGenMassFilter * process.Zprime2muAnalysisSequence)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi')

process.HistosFromPAT.leptonsFromDileptons = True
process.ResolutionUsingMC.leptonsFromDileptons = True
process.ResolutionUsingMC.doQoverP = True

process.p *= rec_level_module(process, process.HistosFromPAT,     'Histos', tracks)
process.p *= rec_level_module(process, process.ResolutionUsingMC, 'Resolution', tracks)

def switch_hlt_name(n):
    process.EfficiencyFromMC.triggerDecision.hltResults = cms.InputTag('TriggerResults', '', n)
    process.HLTSingleObjects.summary = cms.InputTag('hltTriggerSummaryAOD', '', n)
    process.HLTSingleObjects.leptons = [cms.InputTag('hltL3MuonCandidates', '', n)]

#switch_hlt_name('REDIGI38X')

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
        ('dy20',   '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/tucker-effres_dy20-8ca75260210b8943d361f4da5b0c0bcc/USER', 20, 120),
        ('dy120',  '/dy120-HLT-384p3-START38_V12/tucker-effres_mydy120-b62a83c345cd135ef96a2f3fe22d5e32/USER', 120, 200),
        ('dy200',  '/dy200-HLT-384p3-START38_V12/tucker-effres_mydy200-b62a83c345cd135ef96a2f3fe22d5e32/USER', 200, 500),
        ('dy500',  '/dy500-HLT-384p3-START38_V12/tucker-effres_mydy500-b62a83c345cd135ef96a2f3fe22d5e32/USER', 500, 800),
        ('dy800',  '/dy800-HLT-384p3-START38_V12/tucker-effres_mydy800-b62a83c345cd135ef96a2f3fe22d5e32/USER', 800, 20000),
        ('zp1000', '', -20000, 20000),
        ('zp1250', '', -20000, 20000),
        ('zp1500', '', -20000, 20000),
        ('zp1750', '', -20000, 20000),
        ]

    just_testing = 'testing' in sys.argv

    if 'dump_files' in sys.argv:
        for sample in samples:
            print "('%s'," % sample[0]
            os.system('dbss ana02 find file where dataset=%s' % sample[1])
        sys.exit(0)

    for name, dataset, lo, hi in samples:
        if 'zp' in name or name in ['dy20']:
            continue
        
        open('crab.cfg', 'wt').write(crab_cfg % locals())

        new_py = open('histos.py').read()
        new_cut = dy_gen_mass_cut % locals()
        new_py += '\nprocess.DYGenMassFilter.cut = "%(new_cut)s"\n' % locals()
        if name == 'dy20':
            new_py += '\nswitch_hlt_name("REDIGI38X")\n'
        open('histos_crab.py', 'wt').write(new_py)

        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg histos_crab.py')
