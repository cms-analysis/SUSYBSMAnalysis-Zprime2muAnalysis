#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff import loose_cut, trigger_match, tight_cut, allDimuons

process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/datamc_zmumu/5222c20b53e3c47b6c8353d464ee954c/pat_42_3_74A.root']

# Define the numerators and denominators, removing cuts from the
# allDimuons maker. "NoX" means remove cut X (i.e. the denominators),
# "NoNo" means remove nothing (i.e. the numerators). This will break
# if loose_, tight_ cut strings are changed...
process.allDimuonsNoPt      = allDimuons.clone(loose_cut = loose_cut.replace(' && pt > 35.', ''))
process.allDimuonsNoPt35Pt5 = allDimuons.clone(loose_cut = loose_cut.replace(' && pt > 35.', ' && pt > 5.'))
process.allDimuonsNoIso     = allDimuons.clone(loose_cut = loose_cut.replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', ''))
process.allDimuonsNoTkHits  = allDimuons.clone(loose_cut = loose_cut.replace(' && innerTrack.hitPattern.numberOfValidTrackerHits > 10', ''))
process.allDimuonsNoDB      = allDimuons.clone(tight_cut = loose_cut.replace(' && abs(dB) < 0.2', ''))
process.allDimuonsNoGlbChi2 = allDimuons.clone(tight_cut = loose_cut.replace(' && globalTrack.normalizedChi2 < 10', ''))
process.allDimuonsNoPxHits  = allDimuons.clone(tight_cut = loose_cut.replace(' && globalTrack.hitPattern.numberOfValidPixelHits >= 1', ''))
process.allDimuonsNoMuHits  = allDimuons.clone(tight_cut = loose_cut.replace(' && globalTrack.hitPattern.numberOfValidMuonHits > 0', ''))
process.allDimuonsNoMuSegs  = allDimuons.clone(tight_cut = loose_cut.replace(' && numberOfMatchedStations > 1', ''))
process.allDimuonsNoTkMuon  = allDimuons.clone(tight_cut = loose_cut.replace(' && isTrackerMuon', ''))
process.allDimuonsNoTrgMtch = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))
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
process.dimuonsNoB2B     = process.dimuons.clone()
process.dimuonsNoVtxProb = process.dimuons.clone()
delattr(process.dimuonsNoB2B,     'back_to_back_cos_angle_min')
delattr(process.dimuonsNoVtxProb, 'vertex_chi2_max')
process.p *= process.allDimuons
for dimu in ['dimuonsNoB2B', 'dimuonsNoVtxProb']:
    hists = HistosFromPAT.clone(dilepton_src = dimu, leptonsFromDileptons = True)
    setattr(process, dimu.replace('dimuons', ''), hists)
    process.p *= getattr(process, dimu) * hists

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
job_control

[USER]
ui_working_dir = crab/crab_ana_nminus1_%(name)s
return_data = 1
'''

    just_testing = 'testing' in sys.argv
    if not 'no_data' in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2011AMuonsOnly
        Run2011AMuonsOnly.writeJSON('tmp.json')

        dataset_details = [
            ('SingleMu2011A_May10',                '/SingleMu/tucker-datamc_SingleMuRun2011A_May10_new-b2cd34b4395e3cc0cd295229bc3495ca/USER'),
            ('SingleMu2011A_Prompt_165071_165999', '/SingleMu/tucker-datamc_SingleMu2011A_prompt_165071_165999_20110601165210-8788f1b70631d1fb57e97a89f5e8007c/USER'),
            ('SingleMu2011A_Prompt_166000_166562', '/SingleMu/tucker-datamc_SingleMu2011A_prompt_166000_166562_20110608180104-8788f1b70631d1fb57e97a89f5e8007c/USER'),
            ]

        for name, ana_dataset in dataset_details:
            print name

            new_py = open('nminus1effs.py').read()
            new_py += "\nprocess.GlobalTag.globaltag = 'GR_R_42_V13::All'\n"
            open('nminus1effs_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            job_control = '''
total_number_of_lumis = -1
number_of_jobs = 20
lumi_mask = tmp.json'''
            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crab.cfg', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab -create -submit all')

        if not just_testing:
            os.system('rm crab.cfg nminus1effs_crab.py nminus1effs_crab.pyc tmp.json')

    if not 'no_mc' in sys.argv:
        crab_cfg = crab_cfg.replace('job_control','''
total_number_of_events = -1
events_per_job = 50000
''')

        dataset_details = [
            ('zmumu', '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/tucker-datamc_zmumu-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('ttbar', '/TT_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_ttbar-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('dy200', '/DYToMuMu_M-200_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_dy200-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('dy500', '/DYToMuMu_M-500_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_dy500-5222c20b53e3c47b6c8353d464ee954c/USER'),
            #('zssm750', ''),
            ]

        for name, ana_dataset in dataset_details:
            if 'test' in name:
                continue
                crab_cfg = crab_cfg.replace('total_number_of_events = -1', 'total_number_of_events = 55000')

            print name
            open('crab.cfg', 'wt').write(crab_cfg % locals())
            if not just_testing:
                os.system('crab -create -submit all')

        if not just_testing:
            os.system('rm crab.cfg')
