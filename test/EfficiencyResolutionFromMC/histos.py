#!/usr/bin/env python

# Aggggg... too many flags. This workflow should be restructured to
# run like DataMCSpectraComparison/histos.py, where all modes get run
# in the same job.
use_old_selection = False
restrict_mass_window = False
# intime_bin numbering: bin 0 = 0-5, bin 1 = 6-11, bin 2 = 12-26
# late_bin numbering: bin 0 = 0-9, bin 2 = 10-26
intime_bin, late_bin = -1, -1
check_prescaled_path = False
#check_prescaled_path = True
acc_both_24 = False

################################################################################

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.maxEvents.input = -1
#process.source.fileNames = ['file:./pat.root']
process.source.fileNames = [ '/store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/effres_dy2300to3500/160518_094214/0000/pat_10.root',
                            '/store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/effres_dy2300to3500/160518_094214/0000/pat_11.root',
                            '/store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/effres_dy2300to3500/160518_094214/0000/pat_12.root',
                            '/store/user/rradogna/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/effres_dy2300to3500/160518_094214/0000/pat_13.root',
                            ]
process.options.wantSummary = True

ex = ''

if use_old_selection:
    ex += 'oldsel'
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_to_old_selection
    switch_to_old_selection(process)

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
tracks = ['global', 'inner', 'tpfms', 'picky', 'tunep', 'tmr', 'tunepnew']
#tracks = ['tunepnew']
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
process.VBTFEfficiencyFromMC = process.EfficiencyFromMC.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta_2 = 2.4)
process.VBTFEfficiencyFromMCnoTrigger = process.EfficiencyFromMCnoTrigger.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta_2 = 2.4) ### NO TRIGGER PROCESS

# Temporarily disable explicit checks on Level-1 decision until we
# figure out which branch to use.
#for eff in [process.EfficiencyFromMC, process.VBTFEfficiencyFromMC]:
### NO TRIGGER PROCESSES
for eff in [process.EfficiencyFromMCnoTrigger, process.VBTFEfficiencyFromMCnoTrigger]:
    eff.check_l1 = False

if acc_both_24:
    #for eff in [process.EfficiencyFromMC, process.VBTFEfficiencyFromMC]:
    ### NO TRIGGER PROCESSES
    for eff in [process.EfficiencyFromMCnoTrigger, process.VBTFEfficiencyFromMCnoTrigger]:
        eff.acceptance_max_eta_1 = 2.4
    from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match
    for d in [process.allDimuonsVBTF, process.allDimuons]:
        d.tight_cut = trigger_match.replace(' && abs(userFloat("TriggerMatchEta")) < 2.1', '')
    ex += 'accboth24'

#p2 = process.HardInteractionFilter * process.Zprime2muAnalysisSequencePlain * process.EfficiencyFromMC * process.allDimuonsVBTF * process.dimuonsVBTF * process.VBTFEfficiencyFromMC
### NO TRIGGER p2
p2 = process.HardInteractionFilter * process.Zprime2muAnalysisSequencePlain * process.EfficiencyFromMCnoTrigger * process.allDimuonsVBTF * process.dimuonsVBTF * process.VBTFEfficiencyFromMCnoTrigger
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
    #- Calculate efficiency for the thresholds used to count Z events
    from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import prescaled_trigger_pt_threshold, prescaled_offline_pt_threshold
    min_hlt_pt     = prescaled_trigger_pt_threshold
    min_offline_pt = prescaled_offline_pt_threshold
    ex += 'mu' + str(min_hlt_pt)

    l1,hlt = 'L1_SingleMu16er', 'HLT_Mu27_v2'
    from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match

    #for eff in [process.EfficiencyFromMC, process.VBTFEfficiencyFromMC]:
    ### NO TRIGGER
    for eff in [process.EfficiencyFromMCnoTrigger, process.VBTFEfficiencyFromMCnoTrigger]:
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

# Don't care about resolution/etc. when checking the prescaled path.
if not check_prescaled_path:
    process.p  = cms.Path(p)

    # Duplicate all the lepton+dimuon+histogram production for each of
    # what we call "rec levels", i.e. the different TeV reconstructors.
    process.p *= rec_level_module(process, process.HistosFromPAT,     'Histos',     tracks)
    process.p *= rec_level_module(process, process.ResolutionUsingMC, 'Resolution', tracks)

process.p2 = cms.Path(p2)

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()
###############################################################################
    
import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_effres_%(ex)s%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'
#config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20000
#config.Data.unitsPerJob  = 100
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/alfloren'
config.Data.ignoreLocality = True

config.Site.storageSite = 'T2_CH_CERN'
                          
'''
        
    samples = [
               ('dy120','/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/alfloren-effres_dy50to120-2b07f03a08b1e3f01e70977ff823c3f4/USER',50, 120),
               ('dy200','/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/alfloren-effres_dy120to200-2b07f03a08b1e3f01e70977ff823c3f4/USER',120, 200),
               ('dy400','/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/alfloren-effres_dy200to400-2b07f03a08b1e3f01e70977ff823c3f4/USER',200, 400),
               #('dy800','/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/rradogna-effres_dy400to800_s-4acec3425cf1ea1237c0a164f83c7ac2/USER',400, 800),
               ('dy1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/alfloren-effres_dy800to1400-2b07f03a08b1e3f01e70977ff823c3f4/USER',800, 1400),
               ('dy2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/alfloren-effres_dy1400to2300-2b07f03a08b1e3f01e70977ff823c3f4/USER',1400, 2300),
               ('dy3500','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/alfloren-effres_dy1400to2300-2b07f03a08b1e3f01e70977ff823c3f4/USER',2300, 3500),
               ('dy4500','/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/alfloren-effres_dy3500to4500-2b07f03a08b1e3f01e70977ff823c3f4/USER',3500, 4500),
               ('dy6000','/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/alfloren-effres_dy4500to6000-2b07f03a08b1e3f01e70977ff823c3f4/USER',4500, 6000),
               
        ]

    resolutions = {
        500 : 0.048,
        750 : 0.060,
        1000: 0.071,
        1250: 0.080,
        1500: 0.086,
        1750: 0.089,
        2000: 0.089,
        2250: 0.2,
        2500: 0.2,
        2750: 0.2,
        3000: 0.2,
        5000: 0.2
        }
    
    just_testing = 'testing' in sys.argv

    if not restrict_mass_window:
        ex += 'nomasswin'

    if ex:
        ex += '_'
    
    for name, dataset, lo, hi in samples:
        open('crabConfig.py', 'wt').write(crab_cfg % locals())

        if restrict_mass_window and ('zp' in name or 'rs' in name):
            mass = name.replace('zp', '').replace('rs', '')
            mass = mass.replace('_c1', '').replace('_c2', '')
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
            os.system('crab submit  -c crabConfig.py')
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
