#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process
process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-120_7TeV-pythia6/effres_dy120/b62a83c345cd135ef96a2f3fe22d5e32/pat_3_1_Ped.root']
process.options.wantSummary = True

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
tracks = ['global', 'inner', 'tpfms', 'picky', 'pmc', 'tmr', 'sigmaswitch']
rec_levels(process, tracks)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HardInteractionFilter_cfi')
process.HardInteractionFilterRes = process.HardInteractionFilter.clone(use_resonance_mass=True)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')

#process.HardInteractionFilter.use_resonance_mass = True
#process.EfficiencyFromMC.use_resonance_mass_denom = True

process.HLTSingleObjects = cms.EDProducer('HLTLeptonsFromTriggerEvent',
                                          summary = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
                                          leptons = cms.VInputTag(cms.InputTag('hltL3MuonCandidates', '', 'HLT'))
                                          )
process.EfficiencyFromMC.hlt_obj_src = 'HLTSingleObjects'
process.EfficiencyFromMC.hlt_single_min_pt = 15

# Since LooseTightPairSelector ignores the cutFor that
# Zprime2muLeptonProducer sets, don't need to redo leptons for the
# VBTF path.
import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
process.allDimuonsVBTF = VBTFSelection.allDimuons.clone()
process.dimuonsVBTF = VBTFSelection.dimuons.clone(src = 'allDimuonsVBTF')
process.VBTFEfficiencyFromMC = process.EfficiencyFromMC.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta = 2.1)

process.p2 = cms.Path(process.HardInteractionFilter * process.Zprime2muAnalysisSequencePlain * process.HLTSingleObjects * process.EfficiencyFromMC * process.allDimuonsVBTF * process.dimuonsVBTF * process.VBTFEfficiencyFromMC)
process.p = cms.Path(process.HardInteractionFilterRes * process.Zprime2muAnalysisSequence)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi')

process.HistosFromPAT.leptonsFromDileptons = True
process.ResolutionUsingMC.leptonsFromDileptons = True
process.ResolutionUsingMC.doQoverP = True

process.p *= rec_level_module(process, process.HistosFromPAT,     'Histos', tracks)
process.p *= rec_level_module(process, process.ResolutionUsingMC, 'Resolution', tracks)

def switch_hlt_name(n):
    process.EfficiencyFromMC.triggerDecision.hltResults = cms.InputTag('TriggerResults', '', n)
    process.VBTFEfficiencyFromMC.triggerDecision.hltResults = cms.InputTag('TriggerResults', '', n)
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
events_per_job = 60000

[USER]
ui_working_dir = crab/crab_ana_effres_%(name)s
return_data = 1
'''

    samples = [
        ('dy20',  '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/tucker-effres_dy20-8ca75260210b8943d361f4da5b0c0bcc/USER',  20,  60),
        ('dy60',  '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/tucker-effres_dy20-8ca75260210b8943d361f4da5b0c0bcc/USER',  60, 120),
        ('dy120', '/DYToMuMu_M-120_7TeV-pythia6/tucker-effres_dy120-b62a83c345cd135ef96a2f3fe22d5e32/USER',      120, 200),
        ('dy200', '/DYToMuMu_M-200_7TeV-pythia6/tucker-effres_dy200-b62a83c345cd135ef96a2f3fe22d5e32/USER',      200, 500),
        ('dy500', '/DYToMuMu_M-500_7TeV-pythia6/tucker-effres_dy500-b62a83c345cd135ef96a2f3fe22d5e32/USER',      500, 800),
        ('dy800', '/DYToMuMu_M-800_7TeV-pythia6/tucker-effres_dy800-b62a83c345cd135ef96a2f3fe22d5e32/USER',      800, 20000),
        ('zp500',  '/ZprimeSSMToMuMu_M-500_7TeV-pythia6/tucker-effres_zp500-b62a83c345cd135ef96a2f3fe22d5e32/USER',   -20000, 20000),
        ('zp750',  '/ZprimeSSMToMuMu_M-750_7TeV-pythia6/tucker-effres_zp750-b62a83c345cd135ef96a2f3fe22d5e32/USER',   -20000, 20000),
        ('zp1000', '/ZprimeSSMToMuMu_M-1000_7TeV-pythia6/tucker-effres_zp1000-b62a83c345cd135ef96a2f3fe22d5e32/USER', -20000, 20000),
        ('zp1250', '/ZprimeSSMToMuMu_M-1250_7TeV-pythia6/tucker-effres_zp1250-b62a83c345cd135ef96a2f3fe22d5e32/USER', -20000, 20000),
        ('zp1500', '/ZprimeSSMToMuMu_M-1500_7TeV-pythia6/tucker-effres_zp1500-b62a83c345cd135ef96a2f3fe22d5e32/USER', -20000, 20000),
        ('zp1750', '/ZprimeSSMToMuMu_M-1750_7TeV-pythia6/tucker-effres_zp1750-b62a83c345cd135ef96a2f3fe22d5e32/USER', -20000, 20000),
        ]

    resolutions = {
        'zp500' : 0.048,
        'zp750' : 0.060,
        'zp1000': 0.071,
        'zp1250': 0.080,
        'zp1500': 0.086,
        'zp1750': 0.089,
        }
    
    just_testing = 'testing' in sys.argv

    for name, dataset, lo, hi in samples:
        open('crab.cfg', 'wt').write(crab_cfg % locals())

        if False and 'zp' in name:
            mass = float(name.replace('zp', ''))
            res = resolutions[name]
            lo = mass - 1.5*res*mass
            hi = mass + 1.5*res*mass

        new_py = open('histos.py').read()
        new_py += '\nprocess.HardInteractionFilter.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilter.max_mass = "%i"\n' % hi
        new_py += '\nprocess.HardInteractionFilterRes.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilterRes.max_mass = "%i"\n' % hi
        if dataset.startswith('/DYToMuMu_M-20_TuneZ2_7TeV-pythia6'):
            new_py += '\nswitch_hlt_name("REDIGI38X")\n'
        open('histos_crab.py', 'wt').write(new_py)

        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg histos_crab.py histos_crab.pyc')
