#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =[#'file:./pat.root'
#'/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/50000/00105BAD-63E0-E711-8640-02163E0146C5.root',
#'/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F4380BF8-CBCF-E611-8891-0CC47A546E5E.root'
#'/store/mc/RunIIFall17MiniAODv2/ZToEE_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/8A226AAF-AC43-E811-AEF0-0CC47A4D764C.root',
#'file:pickevents.root',
#'/store/data/Run2017E/SingleMuon/MINIAOD/17Nov2017-v1/50000/000DCB8B-2ADD-E711-9100-008CFAF35AC0.root',
#'/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M1300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/DCA196D6-1B79-E711-BADA-90B11C282313.root'
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FEE03939-E778-E711-976D-008CFAF5550C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4E0B4D52-E678-E711-BFE6-0242AC1C0502.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/162AE54B-E678-E711-AFDF-48FD8EE73AC5.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CC2887A8-E778-E711-B955-3417EBE34C27.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F049644B-E678-E711-9036-0025905D1E08.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/20ED2453-E678-E711-9711-0CC47A4C8F30.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/84357A56-E678-E711-8F3E-0CC47A4D764C.root",
#'/store/mc/PhaseIFall16MiniAOD/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/FlatPU28to62HcalNZSRAW_PhaseIFall16_exo52_90X_upgrade2017_realistic_v6_C1-v1/120000/304E419F-CC13-E711-93E9-FA163E0231A1.root',
'/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/70000/12BD4CC4-0751-E811-BCA9-0090FAA58D84.root',
#'/store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam24TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/90F99FA7-E76B-E811-BC35-A0369FD0B362.root'

]                          
process.maxEvents.input = -1
isMC = %(isMC)s
addNTuples = %(addNTuples)s
year = %(year)d
process.GlobalTag.globaltag = '%(GT)s'
process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold, trigger_filters, trigger_path_names, prescaled_trigger_filters, prescaled_trigger_path_names, prescaled_trigger_match_2018, trigger_match_2018

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)
