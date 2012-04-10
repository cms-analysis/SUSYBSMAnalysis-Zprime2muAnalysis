#!/usr/bin/env python

# Aggggg... too many flags. This workflow should be restructured to
# run like DataMCSpectraComparison/histos.py, where all modes get run
# in the same job.
use_old_selection = False
restrict_mass_window = True
# intime_bin numbering: bin 0 = 0-5, bin 1 = 6-11, bin 2 = 12-26
# late_bin numbering: bin 0 = 0-9, bin 2 = 10-26
intime_bin, late_bin = -1, -1 
check_prescaled_path = False
acc_both_24 = True

################################################################################

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.maxEvents.input = 5000
process.source.fileNames = ['/store/user/tucker/ZprimeSSMToMuMu_M-2250_TuneZ2_7TeV-pythia6/effres_zp2250/dd2126535e23ba03e5a28af2e68bf29c/pat_1_1_nHf.root']
process.options.wantSummary = True

ex = ''

if use_old_selection:
    ex += 'oldsel'
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_to_old_selection
    switch_to_old_selection(process)

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
tracks = ['global', 'inner', 'tpfms', 'picky', 'tunep', 'tmr', 'sigmaswitch']
rec_levels(process, tracks)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HardInteractionFilter_cfi')
process.HardInteractionFilterRes = process.HardInteractionFilter.clone(use_resonance_mass=True)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')

#process.HardInteractionFilter.use_resonance_mass = True
#process.EfficiencyFromMC.use_resonance_mass_denom = True

# Since LooseTightPairSelector ignores the cutFor that
# Zprime2muLeptonProducer sets, don't need to redo leptons for the
# VBTF path.
import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
process.allDimuonsVBTF = VBTFSelection.allDimuons.clone()
process.dimuonsVBTF = VBTFSelection.dimuons.clone(src = 'allDimuonsVBTF')
process.VBTFEfficiencyFromMC = process.EfficiencyFromMC.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta_2 = 2.1)

if acc_both_24:
    for eff in [process.EfficiencyFromMC, process.VBTFEfficiencyFromMC]:
        eff.acceptance_max_eta_1 = 2.4
    from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match
    for d in [process.allDimuonsVBTF, process.allDimuons]:
        d.tight_cut = trigger_match.replace(' && abs(userFloat("TriggerMatchEta")) < 2.1', '')
    ex += 'accboth24'

p2 = process.HardInteractionFilter * process.Zprime2muAnalysisSequencePlain * process.EfficiencyFromMC * process.allDimuonsVBTF * process.dimuonsVBTF * process.VBTFEfficiencyFromMC
p  = process.HardInteractionFilterRes * process.Zprime2muAnalysisSequence # this will get all the Histostunep, Histospicky, Histosglobal, etc. below.

if intime_bin in range(0,3) and late_bin in range(0,2):
    # Able to filter on the number of (in-time, late) simulated pileup
    # interactions. By default no filtering is done.
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.GenPileupFilter_cfi') 

    min_intime, max_intime = [(0,5), (6,11), (12,26)][intime_bin]
    min_late, max_late = [(0,9), (10, 26)][late_bin]
    ex += 'PU%i%i%i%i' % (min_intime, max_intime, min_late, max_late)
    process.GenPileupFilter.min_intime = min_intime
    process.GenPileupFilter.max_intime = max_intime
    process.GenPileupFilter.min_late = min_late
    process.GenPileupFilter.max_late = max_late

    p2 = process.GenPileupFilter * p2
    p  = process.GenPileupFilter * p

if check_prescaled_path:
    min_hlt_pt, min_offline_pt = 15, 20
    ex += 'mu' + str(min_hlt_pt)

    l1,hlt = 'L1_SingleMu10', 'HLT_Mu15_v2'
    from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match

    for eff in [process.EfficiencyFromMC, process.VBTFEfficiencyFromMC]:
        eff.triggerDecision.l1Paths = [l1]
        eff.triggerDecision.hltPaths = [hlt]
        eff.hlt_single_min_pt = min_hlt_pt
        eff.acceptance_min_pt = min_offline_pt
        eff.checking_prescaled_path = True

    process.leptons.muon_cuts = 'isGlobalMuon && pt > %i' % min_offline_pt  # Overridden in dimuon construction anyway.

    for d in [process.allDimuonsVBTF, process.allDimuons]:
        assert 'pt > 45' in d.loose_cut.value()
        d.loose_cut = d.loose_cut.value().replace('pt > 45', 'pt > %i' % min_offline_pt)
        if acc_both_24:
            prescaled_trigger_match = prescaled_trigger_match.replace(' && abs(userFloat("TriggerMatchEta")) < 2.1', '')
        d.tight_cut = prescaled_trigger_match

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi')

process.HistosFromPAT.leptonsFromDileptons = True
process.ResolutionUsingMC.leptonsFromDileptons = True
process.ResolutionUsingMC.doQoverP = True

# Don't care about resolution/etc. when checking the Mu15 path.
if not check_prescaled_path:
    process.p  = cms.Path(p)

    # Duplicate all the lepton+dimuon+histogram production for each of
    # what we call "rec levels", i.e. the different TeV reconstructors.
    process.p *= rec_level_module(process, process.HistosFromPAT,     'Histos',     tracks)
    process.p *= rec_level_module(process, process.ResolutionUsingMC, 'Resolution', tracks)

process.p2 = cms.Path(p2)

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
events_per_job = 200000

[USER]
ui_working_dir = crab/crab_ana_effres_%(ex)s%(name)s
return_data = 1
'''
        
    samples = [
        ('dy60',   '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/tucker-effres_dy20-5a92c9dde349d191a78d4eb40ec67d3b/USER',              60,   120),
        ('dy120',  '/DYToMuMu_M-120_TuneZ2_7TeV-pythia6-tauola/tucker-effres_dy120-5a92c9dde349d191a78d4eb40ec67d3b/USER',    120,   200),
        ('dy200',  '/DYToMuMu_M-200_TuneZ2_7TeV-pythia6-tauola/tucker-effres_dy200-5a92c9dde349d191a78d4eb40ec67d3b/USER',    200,   500),
        ('dy500',  '/DYToMuMu_M-500_TuneZ2_7TeV-pythia6-tauola/tucker-effres_dy500-5a92c9dde349d191a78d4eb40ec67d3b/USER',    500,   800),
        ('dy800',  '/DYToMuMu_M-800_TuneZ2_7TeV-pythia6-tauola/tucker-effres_dy800-5a92c9dde349d191a78d4eb40ec67d3b/USER',    800,  1000),
        ('dy1000', '/DYToMuMu_M-1000_TuneZ2_7TeV-pythia6-tauola/tucker-effres_dy1000-5a92c9dde349d191a78d4eb40ec67d3b/USER', 1000, 20000),
        ('zp750',  '/ZprimeSSMToMuMu_M-750_TuneZ2_7TeV-pythia6/tucker-effres_zp750-5a92c9dde349d191a78d4eb40ec67d3b/USER',   -1, 20000),
        ('zp1000', '/ZprimeSSMToMuMu_M-1000_TuneZ2_7TeV-pythia6/tucker-effres_zp1000-5a92c9dde349d191a78d4eb40ec67d3b/USER', -1, 20000),
        ('zp1250', '/ZprimeSSMToMuMu_M-1250_TuneZ2_7TeV-pythia6/tucker-effres_zp1250-5a92c9dde349d191a78d4eb40ec67d3b/USER', -1, 20000),
        ('zp1500', '/ZprimeSSMToMuMu_M-1500_TuneZ2_7TeV-pythia6/tucker-effres_zp1500-5a92c9dde349d191a78d4eb40ec67d3b/USER', -1, 20000),
        ('zp1750', '/ZprimeSSMToMuMu_M-1750_TuneZ2_7TeV-pythia6/tucker-effres_zp1750-5a92c9dde349d191a78d4eb40ec67d3b/USER', -1, 20000),
        ]

    resolutions = {
        500 : 0.048,
        750 : 0.060,
        1000: 0.071,
        1250: 0.080,
        1500: 0.086,
        1750: 0.089,
        }
    
    just_testing = 'testing' in sys.argv

    if not restrict_mass_window:
        ex += 'nomasswin'

    if ex:
        ex += '_'
    
    for name, dataset, lo, hi in samples:
        open('crab.cfg', 'wt').write(crab_cfg % locals())

        if restrict_mass_window and ('zp' in name or 'rs' in name):
            mass = name.replace('zp', '').replace('rs', '')
            mass = int(mass)
            res = resolutions[mass]
            lo = mass - 1.5*res*mass
            hi = mass + 1.5*res*mass

        new_py = open('histos.py').read()
        new_py += '\nprocess.HardInteractionFilter.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilter.max_mass = "%i"\n' % hi
        new_py += '\nprocess.HardInteractionFilterRes.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilterRes.max_mass = "%i"\n' % hi
        open('histos_crab.py', 'wt').write(new_py)

        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm crab.cfg histos_crab.py histos_crab.pyc')
