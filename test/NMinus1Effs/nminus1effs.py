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
# try to check those with a simple string test below.

cuts = [
    ('Pt',      'pt > 45'),
    ('DB',      'abs(dB) < 0.2'),
    ('Iso',     'isolationR03.sumPt / innerTrack.pt < 0.10'),
    ('TkLayers','globalTrack.hitPattern.trackerLayersWithMeasurement > 8'),
    ('PxHits',  'globalTrack.hitPattern.numberOfValidPixelHits >= 1'),
    ('MuHits',  'globalTrack.hitPattern.numberOfValidMuonHits > 0'),
    ('MuMatch', ('numberOfMatchedStations > 1', 'isTrackerMuon')),
    ]

for name, cut in cuts:
    if type(cut) != tuple:
        cut = (cut,)

    lc = loose_cut
    for c in cut:
        if c not in lc:
            raise ValueError('cut "%s" not in cut string "%s"' % (c, lc))
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

# Special case to remove |dB| and B2B cuts simultaneously, as they can
# be correlated (anti-cosmics).
process.allDimuonsNoCosm = process.allDimuons.clone(loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''))
process.dimuonsNoCosm = process.dimuons.clone(src = 'allDimuonsNoCosm')
delattr(process.dimuonsNoCosm, 'back_to_back_cos_angle_min')
process.NoCosm = HistosFromPAT.clone(dilepton_src = 'dimuonsNoCosm', leptonsFromDileptons = True)
process.p *= process.allDimuonsNoCosm * process.dimuonsNoCosm * process.NoCosm

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
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2012MuonsOnly_ll
        Run2012MuonsOnly_ll.writeJSON('tmp.json')

        dataset_details = [
            ('SingleMuRun2012A_13Jul2012_190450_193751', '/SingleMu/slava-datamc_SingleMuRun2012A-13Jul2012_190450_193751_20121011073628-426a2d966f78bce6bde85f3ed41c07ba/USER'),
            ('SingleMuRun2012A_06Aug2012_190782_190949', '/SingleMu/slava-datamc_SingleMuRun2012A-recover-06Aug2012_190782_190949_20121011120430-426a2d966f78bce6bde85f3ed41c07ba/USER'),
            ('SingleMuRun2012B_13Jul2012_193752_196531', '/SingleMu/slava-datamc_SingleMuRun2012B-13Jul2012_193752_196531_20121012044921-426a2d966f78bce6bde85f3ed41c07ba/USER'),
            ('SingleMuRun2012C_24Aug2012_197556_198913', '/SingleMu/slava-datamc_SingleMuRun2012C-24Aug2012_197556_198913_20121012113325-426a2d966f78bce6bde85f3ed41c07ba/USER'),
            ('SingleMuRun2012C_Prompt_198934_203772',    '/SingleMu/slava-datamc_SingleMuRun2012C-Prompt_198934_203772_20121015023300-8627c6a48d2426dec4aa557620a039a0/USER'),
            ('SingleMuRun2012D_Prompt_203773_204563',    '/SingleMu/slava-datamc_SingleMuRun2012D-Prompt_203773_204563_20121016104501-8627c6a48d2426dec4aa557620a039a0/USER'),
            ('SingleMuRun2012D_Prompt_204564_206087',    '/SingleMu/slava-datamc_SingleMuRun2012D-Prompt_204564_206087_20121029121943-8627c6a48d2426dec4aa557620a039a0/USER'),
            ('SingleMuRun2012D-Prompt_206088_206539',    '/SingleMu/slava-datamc_SingleMuRun2012D-Prompt_206088_206539_20121112085341-8627c6a48d2426dec4aa557620a039a0/USER'),
            ]

        for name, ana_dataset in dataset_details:
            print name

            new_py = open('nminus1effs.py').read()
            new_py += "\nprocess.GlobalTag.globaltag = 'GR_P_V42_AN2::All'\n"
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

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
        samples = [zmumu, ttbar, dy120, dy200, dy500, dy800, dy1000, dy1500, dy2000, inclmu15]
        for sample in samples:
            print sample.name
            open('crab.cfg', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab -create -submit all')

        if not just_testing:
            os.system('rm crab.cfg')
