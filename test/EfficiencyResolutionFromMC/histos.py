#!/usr/bin/env python

miniAOD = True

# intime_bin numbering: bin 0 = 0-5, bin 1 = 6-11, bin 2 = 12-26
# late_bin numbering: bin 0 = 0-9, bin 2 = 10-26
intime_bin, late_bin = -1, -1
check_prescaled_path = False


################################################################################

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.maxEvents.input = 1000
process.source.fileNames = [
'/store/mc/PhaseIFall16MiniAOD/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/FlatPU28to62HcalNZSRAW_PhaseIFall16_exo52_90X_upgrade2017_realistic_v6_C1-v1/120000/304E419F-CC13-E711-93E9-FA163E0231A1.root'
                            ]
process.options.wantSummary = True

ex = ''

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import  leptons,leptonsMini, muonPhotonMatchMiniAOD, muonPhotonMatch, allDimuons, dimuons, rec_level_module
#tracks = ['global', 'inner', 'tpfms', 'picky', 'tunep', 'tmr', 'tunepnew']
tracks = ['tunepnew']


process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HardInteractionFilter_cfi')
process.HardInteractionFilterRes = process.HardInteractionFilter.clone(use_resonance_mass=True)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')

import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
process.allDimuonsVBTF = VBTFSelection.allDimuons.clone()
process.dimuonsVBTF = VBTFSelection.dimuons.clone(src = 'allDimuonsVBTF')
process.VBTFEfficiencyFromMCMini = process.EfficiencyFromMCMini.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta_2 = 2.4)
process.VBTFEfficiencyFromMC = process.EfficiencyFromMC.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta_2 = 2.4)
process.VBTFEfficiencyFromMCnoTrigger = process.EfficiencyFromMCnoTrigger.clone(dimuon_src = 'dimuonsVBTF', acceptance_max_eta_2 = 2.4) ### NO TRIGGER PROCESS

# Temporarily disable explicit checks on Level-1 decision until we
# figure out which branch to use.
#for eff in [process.EfficiencyFromMC, process.VBTFEfficiencyFromMC]:
### NO TRIGGER PROCESSES
for eff in [process.EfficiencyFromMCnoTrigger, process.VBTFEfficiencyFromMCnoTrigger, process.VBTFEfficiencyFromMCMini, process.EfficiencyFromMCMini]:
    eff.check_l1 = False


# this will get all the Histostunep, Histospicky, Histosglobal, etc. below.
if miniAOD:
    process.HardInteractionFilterRes.hardInteraction.src = cms.InputTag('prunedGenParticles')
    process.HardInteractionFilter.hardInteraction.src = cms.InputTag('prunedGenParticles')
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
    electrons_miniAOD(process)
    
    
    process.leptons = process.leptonsMini.clone()
    process.Zprime2muAnalysisSequence = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatchMiniAOD * process.egmGsfElectronIDSequence* process.leptons)
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatchMiniAOD * process.egmGsfElectronIDSequence * process.leptons * process.allDimuons * process.dimuons)
   
else:
    process.Zprime2muAnalysisSequence = cms.Sequence(process.HardInteractionFilterRes * process.muonPhotonMatch * process.leptons)
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatch * process.leptons * process.allDimuons * process.dimuons)
   
for t in tracks:
        ad = process.allDimuons.clone()
        label = 'leptons:%s' % t
        setattr(process, 'allDimuons' + t, ad)

        d = process.dimuons.clone()
        d.src = 'allDimuons' + t
        setattr(process, 'dimuons' + t, d)

        process.Zprime2muAnalysisSequence *= ad
        process.Zprime2muAnalysisSequence *= d


p =  process.Zprime2muAnalysisSequence # this will get all the Histostunep, Histospicky, Histosglobal
p2 = process.Zprime2muAnalysisSequencePlain

if miniAOD:
    p2 = p2 * process.EfficiencyFromMCMini * process.allDimuonsVBTF * process.dimuonsVBTF * process.VBTFEfficiencyFromMCMini
#p2 = process.Zprime2muAnalysisSequencePlain * process.EfficiencyFromMCnoTrigger #* process.allDimuonsVBTF make only the full sequence for TuneP and the we appy the efficiency module
else:
    p2 = p2 * process.EfficiencyFromMC * process.allDimuonsVBTF * process.dimuonsVBTF * process.VBTFEfficiencyFromMC

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.HistosFromPAT.leptonsFromDileptons = True
process.HistosFromPAT_MiniAOD.leptonsFromDileptons = True

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi')
if miniAOD:
    process.ResolutionUsingMC.hardInteraction.src = cms.InputTag('prunedGenParticles')
process.ResolutionUsingMC.leptonsFromDileptons = True
process.ResolutionUsingMC.doQoverP = True

if not check_prescaled_path:
    process.p  = cms.Path(p)

    # Duplicate all the lepton+dimuon+histogram production for each of
    # what we call "rec levels", i.e. the different TeV reconstructors.
    #HistosFromPat is applied only on leptons (so no selection applied)
    if miniAOD:
        process.p *= rec_level_module(process, process.HistosFromPAT_MiniAOD,     'Histos',     tracks)
    else:
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
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 200000
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/ferrico'
config.Data.ignoreLocality = True
config.Site.storageSite = 'T2_IT_Bari'
                          
'''
        
    samples = [

	    ('dy50to120', '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_EXO50_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 50., 120.),
    	('dy120to200', '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo51_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 120, 200),
	    ('dy200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo52_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 200, 400),
    	('dy400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo53_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 400, 800),
	    ('dy800to1400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo54_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 800, 1400),
    	('dy1400to2300', '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo55_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 1400, 2300),
	    ('dy2300to3500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo56_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 2300, 3500),
    	('dy3500to4500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo57_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 3500, 4500),
	    ('dy4500to6000', '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/PhaseIFall16MiniAOD-FlatPU28to62HcalNZSRAW_PhaseIFall16_exo58_90X_upgrade2017_realistic_v6_C1-v1/MINIAODSIM', 4500, 6000),
              
               
        ]

    just_testing = 'testing' in sys.argv

    for name, dataset, lo, hi in samples:
        open('crabConfig.py', 'wt').write(crab_cfg % locals())

        new_py = open('histos.py').read()
        new_py += '\nprocess.HardInteractionFilter.min_mass = %i\n' % lo
        new_py += '\nprocess.HardInteractionFilter.max_mass = %i\n' % hi
        new_py += '\nprocess.HardInteractionFilterRes.min_mass = %i\n' % lo
        new_py += '\nprocess.HardInteractionFilterRes.max_mass = %i\n' % hi
        open('histos_crab.py', 'wt').write(new_py)
        
        if not just_testing:
            os.system('crab submit  -c crabConfig.py')
#             os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
