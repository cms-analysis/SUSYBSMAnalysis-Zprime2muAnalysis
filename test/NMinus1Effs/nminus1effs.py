#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection_cff import loose_cut, trigger_match, tight_cut, allDimuons

# Define the numerators and denominators, removing cuts from the
# allDimuons maker. "NoX" means remove cut X (i.e. the denominators),
# "NoNo" means remove nothing (i.e. the numerators). This will break
# if loose_, tight_ cut strings are changed...
process.allDimuonsNoIso     = allDimuons.clone(loose_cut = loose_cut.replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', ''))
process.allDimuonsNoTkHits  = allDimuons.clone(loose_cut = loose_cut.replace(' && innerTrack.hitPattern.numberOfValidTrackerHits >= 10', ''))
process.allDimuonsNoDB      = allDimuons.clone(tight_cut = tight_cut.replace('dB < 0.2 && ', ''))
process.allDimuonsNoGlbChi2 = allDimuons.clone(tight_cut = tight_cut.replace(' && globalTrack.normalizedChi2 < 10', ''))
process.allDimuonsNoPxHits  = allDimuons.clone(tight_cut = tight_cut.replace(' && innerTrack.hitPattern.numberOfValidPixelHits >= 1', ''))
process.allDimuonsNoMuStns  = allDimuons.clone(tight_cut = tight_cut.replace(' && globalTrack.hitPattern.muonStationsWithValidHits >= 2', ''))
process.allDimuonsNoTkMuon  = allDimuons.clone(tight_cut = tight_cut.replace(' && isTrackerMuon', ''))
process.allDimuonsNoTrgMtch = allDimuons.clone(tight_cut = tight_cut.replace(' && ' + trigger_match, ''))
process.allDimuonsNoNo      = allDimuons.clone()

alldimus = [x for x in dir(process) if 'allDimuonsNo' in x]
process.p = cms.Path(process.goodDataFilter * process.muonPhotonMatch * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))

# For all the allDimuons producers, make dimuons producers, and
# analyzers to make the histograms.
for alld in alldimus:
    dimu = process.dimuons.clone(src = alld)
    name = alld.replace('allD', 'd')
    setattr(process, name, dimu)
    hists = HistosFromPAT.clone(dilepton_src = name, leptonsFromDileptons = True)
    setattr(process, name.replace('dimuons', ''), hists)
    process.p *= dimu * hists

# Handle the cuts that have to be applied at the
# Zprime2muCompositeCandidatePicker level.
process.p *= process.allDimuons
process.dimuonsNoB2B     = process.dimuons.clone(back_to_back_cos_angle_min = -1)
process.dimuonsNoVtxProb = process.dimuons.clone(vertex_chi2_max = -1)

for dimu in ['dimuonsNoB2B', 'dimuonsNoVtxProb']:
    hists = HistosFromPAT.clone(dilepton_src = dimu, leptonsFromDileptons = True)
    setattr(process, dimu.replace('dimuons', ''), hists)
    process.p *= getattr(process, dimu) * hists

if 'data' in sys.argv:
    process.source.fileNames = ['file:crab/crab_datamc_Run2010A/res/merged.root', 'file:crab/crab_datamc_promptB_all/res/merged.root']
    process.GlobalTag.globaltag = 'GR10_P_V10::All'
    from goodlumis import Run2010ABMuonsOnly
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*Run2010ABMuonsOnly)
    process.TFileService.fileName = 'ana_nminus1_data.root'

if 'debugcuts' in sys.argv:
    process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-800_7TeV-pythia6/datamc_dy800/20d63349d4f532c6f4f93a8966ef6c34/pat_2_1_rIJ.root']
    process.maxEvents.input = 1000
    process.GlobalTag.globaltag = 'START38_V14::All'
    process.MessageLogger.categories.append('PrintEvent')
    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.p *= process.triggerSummaryAnalyzerAOD
    process.PrintEvent = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('dimuonsNoNo'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.PrintEventNoTrgMtch = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('dimuonsNoTrgMtch'))
    process.p *= process.PrintEvent * process.PrintEventNoTrgMtch
    
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(ana_dataset)s
dbs_url = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
pset = nminus1effs.py
get_edm_output = 1
total_number_of_events = -1
events_per_job = 20000

[USER]
ui_working_dir = crab/crab_ana_nminus1_%(name)s
return_data = 1
'''

    just_testing = 'testing' in sys.argv
    
    x = [
        ('zmumu', '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/tucker-datamc_zmumu-b4341788d83565203f0d6250b5475e6e/USER'),
        ('zssm750', '/ZprimeSSMToMuMu_M-750_7TeV-pythia6/tucker-datamc_zssm750-b4341788d83565203f0d6250b5475e6e/USER'),
        ('ttbar', '/TTJets_TuneZ2_7TeV-madgraph-tauola/tucker-datamc_ttbar-b4341788d83565203f0d6250b5475e6e/USER'),
        ('dy120', '/DYToMuMu_M-120_7TeV-pythia6/tucker-effres_dy120-b62a83c345cd135ef96a2f3fe22d5e32/USER'),
        ]
    
    for name, ana_dataset in x:
        if name != 'dy120':
            continue
        print name
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not just_testing:
            os.system('crab -create -submit all')
        
    if not just_testing:
        os.system('rm crab.cfg')
