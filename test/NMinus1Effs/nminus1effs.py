#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff import loose_cut, trigger_match, tight_cut, allDimuons

process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/datamc_zmumu/5222c20b53e3c47b6c8353d464ee954c/pat_42_3_74A.root']
process.maxEvents.input = 1000

# Define the numerators and denominators, removing cuts from the
# allDimuons maker. "NoX" means remove cut X entirely (i.e. the
# loose_cut denominators), "TiX" means move cut X from the loose_cut
# to the tight_cut (meaning only one muon instead of two has to pass
# the cut).  "NoNo" means remove nothing (i.e. the numerator). This
# will break if loose_, tight_cut strings are changed upstream, so we
# try to check those with a silly hash next. (Don't want to check full
# string equality since then we have one more string to maintain
# besides the N of the below.)

assert hash(loose_cut) == -5604570599357377599
assert hash(tight_cut) == -2883478064365267914

cuts = [
    ('Pt',      'pt > 35.'),
    ('DB',      'abs(dB) < 0.2'),
    ('GlbChi2', 'globalTrack.normalizedChi2 < 10'),
    ('Iso',     'isolationR03.sumPt / innerTrack.pt < 0.10'),
    ('TkHits',  'globalTrack.hitPattern.numberOfValidTrackerHits > 10'),
    ('PxHits',  'globalTrack.hitPattern.numberOfValidPixelHits >= 1'),
    ('MuHits',  'globalTrack.hitPattern.numberOfValidMuonHits > 0'),
    ('MuMatch', ('numberOfMatchedStations > 1', 'isTrackerMuon')),
    ]

for name, cut in cuts:
    if type(cut) != tuple:
        cut = (cut,)

    lc = loose_cut
    for c in cut:
        lc = lc.replace(' && ' + c, '') # Relies on none of the cuts above being first in the list.

    obj_no = allDimuons.clone(loose_cut = lc)
    setattr(process, 'allDimuonsNo' + name, obj_no)
    
    obj_ti = obj_no.clone(tight_cut = tight_cut + ' && ' + ' && '.join(cut))
    setattr(process, 'allDimuonsTi' + name, obj_ti)

process.allDimuonsNoNo      = allDimuons.clone()
process.allDimuonsNoTrgMtch = allDimuons.clone(tight_cut = tight_cut.replace(trigger_match, ''))

alldimus = [x for x in dir(process) if 'allDimuonsNo' in x or 'allDimuonsTi' in x]

# Sanity check that the replaces above did something.
for x in alldimus:
    if 'NoNo' in x:
        continue
    o = getattr(process, x)
    assert o.loose_cut.value() != loose_cut or o.tight_cut.value() != tight_cut

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
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2011AMuonsOnly_ll
        Run2011AMuonsOnly_ll.writeJSON('tmp.json')

        dataset_details = [
            ('SingleMu2011A_May10',                '/SingleMu/tucker-datamc_SingleMuRun2011A_May10_new-b2cd34b4395e3cc0cd295229bc3495ca/USER'),
            ('SingleMu2011A_Prompt_165071_165999', '/SingleMu/tucker-datamc_SingleMu2011A_prompt_165071_165999_20110601165210-8788f1b70631d1fb57e97a89f5e8007c/USER'),
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
            ('zmumu',    '/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/tucker-datamc_zmumu-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('ttbar',    '/TT_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_ttbar-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('dy120',    '/DYToMuMu_M-120_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_dy120-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('dy200',    '/DYToMuMu_M-200_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_dy200-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('dy500',    '/DYToMuMu_M-500_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_dy500-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ('dy1000',   '/DYToMuMu_M-1000_TuneZ2_7TeV-pythia6-tauola/tucker-datamc_dy1000-5222c20b53e3c47b6c8353d464ee954c/USER'),
            ]

        for name, ana_dataset in dataset_details:
            print name
            open('crab.cfg', 'wt').write(crab_cfg % locals())
            if not just_testing:
                os.system('crab -create -submit all')

        if not just_testing:
            os.system('rm crab.cfg')
