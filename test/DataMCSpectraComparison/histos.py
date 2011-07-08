#!/usr/bin/env python

import sys, os, glob, FWCore.ParameterSet.Config as cms
from pprint import pprint
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
process.maxEvents.input = 5000
process.source.fileNames = ['/store/user/tucker/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/datamc_zmumu/5222c20b53e3c47b6c8353d464ee954c/pat_42_3_74A.root']

# This is needed for random numbers when doing the extra prescaling
# for the Mu15 path.
process.RandomNumberGeneratorService = cms.Service('RandomNumberGeneratorService', PrescaleToCommon = cms.PSet(initialSeed = cms.untracked.uint32(1219)))
Mu15s = ['HLT_Mu15_v2', 'HLT_Mu15_v3', 'HLT_Mu15_v4', 'HLT_Mu15_v5', 'HLT_Mu15_v6']

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match
mu15_trigger_match = ' || '.join('!triggerObjectMatchesByPath("%s").empty()' % p for p in Mu15s)

# Since Mu15 comes with different prescales in different runs/lumis,
# this filter prescales it to a common factor to make things simpler.
process.PrescaleToCommon = cms.EDFilter('PrescaleToCommon',
                                        hlt_process_name = cms.string('HLT'),
                                        trigger_paths = cms.vstring(*Mu15s),
                                        overall_prescale = cms.int32(350),
                                        assume_simulation_has_prescale_1 = cms.bool(True) # Current PAT tuples of MC samples don't have both L1 branches :-(
                                        )

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT

from SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi import dy_gen_mass_cut
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi')

import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionOld_cff as OurSelectionOld
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew

# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
common_dil_cut = ''
dils = [
    ('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    ('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
    ('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
    ('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
    ('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
    ('MuonsPlusElectronsMinus',      '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == -2'),
    ('MuonsMinusElectronsPlus',      '%(leptons_name)s:muons@- %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == 2'),
    ('MuonsPlusElectronsPlus',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == -24'),
    ('MuonsMinusElectronsMinus',     '%(leptons_name)s:muons@- %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == 24'),
    ('MuonsElectronsOppSign',        '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     ''),
    ('MuonsElectronsSameSign',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
    ('MuonsElectronsAllSigns',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
    ]

# Define groups of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add hltFilter
# somewhere below.
cuts = {
    'VBTF'     : VBTFSelection,
    'OurOld'   : OurSelectionOld,
    'OurNew'   : OurSelectionNew,
    'OurNoIso' : OurSelectionNew,
    'EmuVeto'  : OurSelectionNew,
    'Simple'   : OurSelectionNew, # the selection cuts in the module listed here are ignored below
    'VBTFMu15' : VBTFSelection,
    'OurMu15'  : OurSelectionNew,
    }

for cut_name, Selection in cuts.iteritems():
    # Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'Mu15' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > 35', 'pt > 20')
    else:
        muon_cuts = Selection.loose_cut
    leptons = process.leptons.clone(muon_cuts = muon_cuts)
    if cut_name == 'EmuVeto':
        leptons.electron_muon_veto_dR = 0.1
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the Mu15 paths, we don't care about e-mu events.
        if 'Mu15' in cut_name and 'Electron' in dil_name:
            continue
        
        name = cut_name + dil_name
        allname = 'all' + name

        if common_dil_cut and dil_cut:
            dil_cut = common_dil_cut + ' && (%s)' % dil_cut
            
        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)
        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        # Implement the differences to the selections.
        if cut_name == 'Simple':
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.do_remove_overlap = False
            delattr(dil, 'back_to_back_cos_angle_min')
            delattr(dil, 'vertex_chi2_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
        elif 'Mu15' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace('pt > 35', 'pt > 20')
            assert alldil.tight_cut == trigger_match
            alldil.tight_cut = mu15_trigger_match

        # Histos now just needs to know which leptons and dileptons to use.
        histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)

    # Finally, make the path for this set of cuts. Don't use hltFilter
    # here, but rely on the selection to take care of it -- easy way
    # to handle the changing trigger names.
    pathname = 'path' + cut_name
    pobj = process.muonPhotonMatch * reduce(lambda x,y: x*y, path_list)
    if 'VBTF' not in cut_name and cut_name != 'Simple':
        pobj = process.goodDataFilter * pobj
    if 'Mu15' in cut_name:
        pobj = process.PrescaleToCommon * pobj
    path = cms.Path(pobj)
    setattr(process, pathname, path)

def ntuplify(process, hlt_process_name='HLT', fill_gen_info=False):
    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler',
                                           hlt_src = cms.InputTag('TriggerResults', '', hlt_process_name),
                                           dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlinePrimaryVertices')
                                           )
    process.SimpleNtuplerEmu = process.SimpleNtupler.clone(dimu_src = cms.InputTag('SimpleMuonsElectronsAllSigns'))

    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        process.SimpleNtupler.hardInteraction = hardInteraction

    process.pathSimple *= process.SimpleNtupler * process.SimpleNtuplerEmu

def printify(process, hlt_process_name='HLT'):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD', '', hlt_process_name)

    process.PrintEvent = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    process.PrintEventSS = process.PrintEvent.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    process.PrintEventEmu = process.PrintEvent.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    process.pathOurNew *= process.PrintEvent * process.PrintEventSS * process.PrintEventEmu

    process.PrintEventVBTF = process.PrintEvent.clone(dilepton_src = cms.InputTag('VBTFMuonsPlusMuonsMinus'))
    process.pathVBTF *= process.PrintEventVBTF

    process.PrintEventSimple = process.PrintEvent.clone(dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    process.pathSimple *= process.triggerSummaryAnalyzerAOD * process.PrintEventSimple

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.CheckPrescale = cms.EDAnalyzer('CheckPrescale',
                                           hlt_process_name = cms.string(hlt_process_name),
                                           trigger_paths = cms.vstring(*trigger_paths)
                                           )
    process.pCheckPrescale = cms.Path(process.CheckPrescale)

if 'gogo' in sys.argv:
    ntuplify(process) #, fill_gen_info=True)
    #printify(process)
    
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(ana_dataset)s
dbs_url = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
pset = histos_crab.py
get_edm_output = 1
job_control

[USER]
ui_working_dir = crab/crab_ana_datamc_%(name)s
return_data = 1
'''

    just_testing = 'testing' in sys.argv
    
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [
            ('SingleMu2011A_May10',                '/SingleMu/tucker-datamc_SingleMu2011A_May10-3c88448713b4112b83eb5e163e8441f1/USER'),
            ('SingleMu2011A_Prompt_165071_167913', '/SingleMu/tucker-datamc_SingleMu2011A_prompt_165071_167913_20110705115318-e0e58cf0dbd55d2562f61b8061f4c446/USER'),
            ]

        lumi_lists = [
            'Run2011A',
            'Run2011AMuonsOnly',
            #'Run2011APlusDCSOnlyMuonsOnly',
            ]

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll')
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = 'lumi_mask = %s' % json_fn
                
            name = '%s_%s' % (lumi_name, dataset_name)
            print name

            new_py = open('histos.py').read()
            new_py += "\nntuplify(process)\n"
            new_py += "\nprocess.GlobalTag.globaltag = 'GR_R_42_V13::All'\n"
            new_py += "\ncheck_prescale(process, ['HLT_Mu30_v1', 'HLT_Mu30_v2', 'HLT_Mu30_v3', 'HLT_Mu30_v4', 'HLT_Mu30_v5'])\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
total_number_of_lumis = -1
lumis_per_job = 500
%(lumi_mask)s''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crab.cfg', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab -create -submit all')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.cfg'
                print cmd
                os.system(cmd)

        if not just_testing:
            os.system('rm crab.cfg histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
total_number_of_events = -1
events_per_job = 100000
    ''')

        import samples
        combine_dy_samples = len([x for x in samples.samples if x.name in ['dy200', 'dy500', 'dy800', 'dy1000']]) > 0
        print 'combine_dy_samples:', combine_dy_samples

        from samples import samples
        for sample in reversed(samples):
            if 'Summer11' not in sample.dataset:
                continue
            print sample.name

            new_py = open('histos.py').read()
            new_py += "\nprocess.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', '%(hlt_process_name)s')\n" % sample

            sample.fill_gen_info = sample.name in ['zmumu', 'dy120', 'dy200', 'dy500', 'dy800', 'dy1000', 'zssm1000']
            new_py += "\nntuplify(process, hlt_process_name='%(hlt_process_name)s', fill_gen_info=%(fill_gen_info)s)\n" % sample

            if combine_dy_samples and (sample.name == 'zmumu' or 'dy' in sample.name):
                mass_limits = {
                    'zmumu' : (  20,    120),
                    'dy120' : ( 120,    200),
                    'dy200' : ( 200,    500),
                    'dy500' : ( 500,    800),
                    'dy800' : ( 800,   1000),
                    'dy1000': (1000, 100000),
                    }
                lo,hi = mass_limits[sample.name]
                new_cut = dy_gen_mass_cut % locals()
                new_py += '\nprocess.DYGenMassFilter.cut = "%(new_cut)s"\n' % locals()
                new_py += '\nfor pn,p in process.paths.items():\n  setattr(process, pn, cms.Path(process.DYGenMassFilter*p._seq))\n'

            open('histos_crab.py', 'wt').write(new_py)

            open('crab.cfg', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab -create -submit all')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.cfg'
                print cmd
                os.system(cmd)

        if not just_testing:
            os.system('rm crab.cfg histos_crab.py histos_crab.pyc')
