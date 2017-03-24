#!/usr/bin/env python

miniAOD = True

# intime_bin numbering: bin 0 = 0-5, bin 1 = 6-11, bin 2 = 12-26
# late_bin numbering: bin 0 = 0-9, bin 2 = 10-26
intime_bin, late_bin = -1, -1
check_prescaled_path = False


################################################################################

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.maxEvents.input = -1
process.source.fileNames = [#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/02244373-7E03-E611-B581-003048F5B2B4.root',
 			'file:/u/user/msoh/Samples/MC/DY/M6000/0A222CC8-16C7-E611-BFC2-FA163E39B5B6.root', #M4500-6000                    
			'file:/u/user/msoh/Samples/MC/DY/M6000/52FEE967-17C7-E611-9A69-FA163E2EED86.root',
			'file:/u/user/msoh/Samples/MC/DY/M6000/A4999F0B-17C7-E611-A246-1866DAEA6D0C.root',
			'file:/u/user/msoh/Samples/MC/DY/M6000/A8CCA36B-17C7-E611-8E4A-0090FAA57D64.root',
			'file:/u/user/msoh/Samples/MC/DY/M6000/B86337F0-16C7-E611-BA6A-24BE05C44B91.root',
			'file:/u/user/msoh/Samples/MC/DY/M6000/EA325C83-17C7-E611-AFA4-001E67E5E8B6.root',
			'file:/u/user/msoh/Samples/MC/DY/M6000/F214122A-17C7-E611-9134-001E674FC800.root',

			#'file:/u/user/msoh/Data/MC/E412A19B-DCCF-E611-8050-0CC47A546E5E.root',  ##Moriond17 M50-120
			#'file:/u/user/msoh/Data/MC/824C363B-0AC8-E611-B4A5-20CF3027A580.root', ##Moriond17
			#'file:/u/user/msoh/Data/MC/52FEE967-17C7-E611-9A69-FA163E2EED86.root', ##Moriond17 M4500-6000

				#'file:/u/user/msoh/data/mc/m50_120/1A6B76DF-153B-E611-BEC5-0CC47A4DEDF8.root',


                            ]
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000


from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import  leptons,leptonsMini, muonPhotonMatchMiniAOD, muonPhotonMatch, allDimuons, dimuons, rec_level_module

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HardInteractionFilter_cfi')
process.HardInteractionFilterRes = process.HardInteractionFilter.clone(use_resonance_mass=True)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection
process.allDimuonsOur = OurSelection.allDimuons.clone()
process.dimuonsOur = OurSelection.dimuons.clone(src = 'allDimuonsOur')
process.OurEfficiencyFromMCMini = process.EfficiencyFromMCMini.clone(dimuon_src = 'dimuonsOur', acceptance_max_eta_2 = 2.4)
process.OurEfficiencyFromMC = process.EfficiencyFromMC.clone(dimuon_src = 'dimuonsOur', acceptance_max_eta_2 = 2.4)
process.OurEfficiencyFromMCnoTrigger = process.EfficiencyFromMCnoTrigger.clone(dimuon_src = 'dimuonsOur', acceptance_max_eta_2 = 2.4) ### NO TRIGGER PROCESS

# Temporarily disable explicit checks on Level-1 decision until we
# figure out which branch to use.
#for eff in [process.EfficiencyFromMC, process.OurEfficiencyFromMC]:
### NO TRIGGER PROCESSES
for eff in [process.EfficiencyFromMCnoTrigger, process.OurEfficiencyFromMCnoTrigger, process.OurEfficiencyFromMCMini, process.EfficiencyFromMCMini]:
    eff.check_l1 = False


# this will get all the Histostunep, Histospicky, Histosglobal, etc. below.
if miniAOD:
    process.HardInteractionFilterRes.hardInteraction.src = cms.InputTag('prunedGenParticles')
    process.HardInteractionFilter.hardInteraction.src = cms.InputTag('prunedGenParticles')
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
    electrons_miniAOD(process)
    
    
    process.leptons = process.leptonsMini.clone()
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatchMiniAOD * process.egmGsfElectronIDSequence * process.leptons * process.allDimuons * process.dimuons)
   
else:
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatch * process.leptons * process.allDimuons * process.dimuons)
   

p2 = process.Zprime2muAnalysisSequencePlain

if miniAOD:
    p2 = p2 * process.EfficiencyFromMCMini * process.allDimuonsOur * process.dimuonsOur * process.OurEfficiencyFromMCMini

else:
    p2 = p2 * process.EfficiencyFromMC * process.allDimuonsOur * process.dimuonsOur * process.OurEfficiencyFromMC

process.p2 = cms.Path(p2)

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()

###############################################################################

ex = '0210'
    
import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'eff_%(ex)s_%(name)s'
config.General.workArea = 'crab'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_eff.py'
#config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20000
#config.Data.unitsPerJob  = 100
config.Data.publication = False
config.Data.outputDatasetTag = 'eff_%(ex)s_%(name)s'
config.Data.outLFNDirBase = '/store/user/moh'
config.Data.ignoreLocality = True

config.Site.storageSite = 'T2_KR_KNU'
                          
'''
        
    samples = [
             #  ('dy120','/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',50, 120),
             #  ('dy200','/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',120, 200),
             #  ('dy400','/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',200, 400),
             #  ('dy800','/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',400, 800),
             #  ('dy1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',800, 1400),
             #  ('dy2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',1400, 2300),
             #  ('dy3500','/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',2300,3500),
             #  ('dy4500','/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',3500,4500),
             #  ('dy6000','/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',4500,6000),
             #  ('dyInf','/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',6000,-1),
     
                ('dy120','/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',50,120),
	       	#('dy200','/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',120,200),
                #('dy400','/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',200,400),
                #('dy800','/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',400,800),
                #('dy1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',800,1400),
                #('dy2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',1400,2300),
                #('dy3500','/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',2300,3500),
                #('dy4500','/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3500,4500),
                #('dy6000','/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',4500,6000),
                #('dyInf','/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6000,100000),

	   ]

    just_testing = 'testing' in sys.argv

    for name, dataset, lo, hi in samples:
        open('crabConfig.py', 'wt').write(crab_cfg % locals())

        new_py = open('histos.py').read()
        new_py += '\nprocess.HardInteractionFilter.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilter.max_mass = "%i"\n' % hi
        new_py += '\nprocess.HardInteractionFilterRes.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilterRes.max_mass = "%i"\n' % hi
        open('histos_crab.py', 'wt').write(new_py)
        
        if not just_testing:
            os.system('crab submit  -c crabConfig.py')
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
