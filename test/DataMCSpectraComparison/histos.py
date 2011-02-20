#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
process.source.fileNames = ['/store/user/tucker/TTJets_TuneZ2_7TeV-madgraph-tauola/datamc_ttbar/b4341788d83565203f0d6250b5475e6e/pat_9_1_0F5.root']

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT

from SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi import dy_gen_mass_cut
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi')

import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection_cff as OurSelection
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
    'VBTF'  : VBTFSelection,
    'Our'   : OurSelection,
    'OurNew': OurSelectionNew,
    'Simple': OurSelection,
    }

for cut_name in cuts.keys():
    # Keep track of modules to put in the path for this set of cuts.
    path_list = []

    Selection = cuts[cut_name]

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.
    leptons_name = cut_name + 'Leptons'
    leptons = process.leptons.clone(muon_cuts = '' if cut_name == 'Simple' else Selection.loose_cut)
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        name = cut_name + dil_name
        allname = 'all' + name

        if common_dil_cut and dil_cut:
            dil_cut = common_dil_cut + ' && (%s)' % dil_cut
            
        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)
        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        if cut_name == 'Simple':
            alldil.loose_cut = 'isGlobalMuon && (pt > 20. || innerTrack.pt > 20.)'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.do_remove_overlap = False
            delattr(dil, 'back_to_back_cos_angle_min')
            delattr(dil, 'vertex_chi2_max')
                    
        histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

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
    path = cms.Path(pobj)
    setattr(process, pathname, path)

def ntuplify(process, hlt_process_name='HLT'):
    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler', hlt_src = cms.InputTag('TriggerResults', '', hlt_process_name), dimu_src = cms.InputTag('SimpleMuonsAllSigns'))
    process.SimpleNtuplerEmu = cms.EDAnalyzer('SimpleNtupler', hlt_src = cms.InputTag('TriggerResults', '', hlt_process_name), dimu_src = cms.InputTag('SimpleMuonsElectronsAllSigns'))
    process.pathSimple *= process.SimpleNtupler * process.SimpleNtuplerEmu

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')
    process.PrintEvent = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('OurMuonsPlusMuonsMinus'))
    process.PrintEventSS = process.PrintEvent.clone(dilepton_src = cms.InputTag('OurMuonsSameSign'))
    process.PrintEventEmu = process.PrintEvent.clone(dilepton_src = cms.InputTag('OurMuonsElectronsOppSign'))
    process.pathOur *= process.PrintEvent * process.PrintEventSS * process.PrintEventEmu

    process.PrintEventVBTF = process.PrintEvent.clone(dilepton_src = cms.InputTag('VBTFMuonsPlusMuonsMinus'))
    process.pathVBTF *= process.PrintEventVBTF

    process.PrintEventSimple = process.PrintEvent.clone(dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    process.pathSimple *= process.PrintEventSimple

if 'bigfatbean' in sys.argv:
    ntuplify(process)
    for c in ['VBTF', 'Our', 'OurNew']:
        delattr(process, 'path%s' % c)
    process.pathSimple=cms.Path(process.muonPhotonMatch*process.SimpleLeptons*process.allSimpleMuonsAllSigns*process.SimpleMuonsAllSigns*process.SimpleMuonsAllSignsHistos*process.SimpleNtupler)
    
if 'data' in sys.argv:
    process.source.fileNames = ['file:crab/crab_datamc_Run2010A/res/merged.root', 'file:crab/crab_datamc_promptB_all/res/merged.root']
    process.TFileService.fileName = 'ana_datamc_data.root'
    process.GlobalTag.globaltag = 'GR10_P_V10::All'

    from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Run2010AB, Run2010ABMuonsOnly
    if 'all_good' in sys.argv:
        process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*Run2010AB)
    else:
        process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*Run2010ABMuonsOnly)

    printify(process)
    ntuplify(process)

if 'nov4' in sys.argv:
    process.source.fileNames = ['file:crab/crab_datamc_Run2010A_DileptonMu/res/merged.root', 'file:crab/crab_datamc_Run2010B_DileptonMu/res/merged.root']
    process.TFileService.fileName = 'ana_datamc_data.root'
    process.GlobalTag.globaltag = 'GR_R_38X_V15::All'

    from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Nov4Run2010AB, Nov4Run2010ABMuonsOnly
    if 'all_good' in sys.argv:
        process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*Nov4Run2010AB)
    else:
        process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*Nov4Run2010ABMuonsOnly)

    printify(process)
    ntuplify(process)

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
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import Nov4Run2010AB_ll, Nov4Run2010ABMuonsOnly_ll
        x = [
            ('Run2010AAllGood',   '/Mu/tucker-datamc_Run2010A_DileptonMu-4b90be408f306fdc739fecf72d09b336/USER', Nov4Run2010AB_ll),
            ('Run2010BAllGood',   '/Mu/tucker-datamc_Run2010B_DileptonMu-4b90be408f306fdc739fecf72d09b336/USER', Nov4Run2010AB_ll),
            ('Run2010AMuonsOnly', '/Mu/tucker-datamc_Run2010A_DileptonMu-4b90be408f306fdc739fecf72d09b336/USER', Nov4Run2010ABMuonsOnly_ll),
            ('Run2010BMuonsOnly', '/Mu/tucker-datamc_Run2010B_DileptonMu-4b90be408f306fdc739fecf72d09b336/USER', Nov4Run2010ABMuonsOnly_ll),
            ]

        for name, ana_dataset, lumi_list in x:
            print name
            new_py = open('histos.py').read()
            new_py += "\nntuplify(process)\n"
            new_py += "\nprocess.GlobalTag.globaltag = 'GR_R_38X_V15::All'\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()
            lumi_list.writeJSON('tmp.json')

            new_crab_cfg = new_crab_cfg.replace('job_control','''
total_number_of_lumis = -1
number_of_jobs = 5
lumi_mask = tmp.json''')

            open('crab.cfg', 'wt').write(new_crab_cfg)
            if not just_testing:
                os.system('crab -create -submit all')

        if not just_testing:
            os.system('rm crab.cfg histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
total_number_of_events = -1
events_per_job = 20000
    ''')

        from samples import samples
        for sample in samples:
            print sample.name

            new_py = open('histos.py').read()
            new_py += "\nprocess.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', '%(hlt_process_name)s')\n" % sample
            new_py += "\nntuplify(process, hlt_process_name='%(hlt_process_name)s')\n" % sample

            if sample.name == 'zmumu' or 'dy' in sample.name:
                mass_limits = {
                    'zmumu': ( 20, 200),
                    'dy200': (200, 500),
                    'dy500': (500, 800),
                    'dy800': (800, 100000),
                    }
                lo,hi = mass_limits[sample.name]
                new_cut = dy_gen_mass_cut % locals()
                new_py += '\nprocess.DYGenMassFilter.cut = "%(new_cut)s"\n' % locals()
                new_py += '\nfor pn,p in process.paths.items():\n  setattr(process, pn, cms.Path(process.DYGenMassFilter*p._seq))\n'

            open('histos_crab.py', 'wt').write(new_py)

            open('crab.cfg', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab -create -submit all')

        if not just_testing:
            os.system('rm crab.cfg histos_crab.py histos_crab.pyc')
