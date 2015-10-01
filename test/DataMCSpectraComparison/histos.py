#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

process.source.fileNames =['file:./pat.root',
                           
#                           '/store/user/rradogna/SingleMuon/datamc_SingleMuonRun2015C-Prompt_253888_254914_20150831150018/150831_130042/0000/pat_1.root',
#

                           ]
process.maxEvents.input = -1
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'## solo per proare i dati
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A'## solo per proare i mc
#process.GlobalTag.globaltag = '74X_mcRun2_startup_realistic50ns_v0' #mc startup
#process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
#HistosFromPAT.useMadgraphWeight= cms.untracked.bool(False) #if we want to re-wigth madgraph samples. For other samples the weight is always 1. Does not work with data.

# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
#import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionOld_cff as OurSelectionOld
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2011EPS_cff as OurSelection2011EPS
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff as OurSelectionDec2012

# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
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

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
#    'VBTF'     : VBTFSelection,
#    'OurOld'   : OurSelectionOld,
#    'OurEPS'   : OurSelection2011EPS,
    #'OurNew'   : OurSelectionNew,
    'Our2012'  : OurSelectionDec2012,
    #'OurNoIso' : OurSelectionDec2012,
    #'EmuVeto'  : OurSelectionDec2012,
    'Simple'   : OurSelectionDec2012, # The selection cuts in the module listed here are ignored below.
#    'VBTFMuPrescaled' : VBTFSelection,
    #'OurMuPrescaledNew'  : OurSelectionNew,
#    'OurMuPrescaled2012' : OurSelectionDec2012
    }

# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
    # Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present.
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'MuPrescaled' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut
    leptons = process.leptons.clone(muon_cuts = muon_cuts)
    if cut_name == 'EmuVeto':
        leptons.electron_muon_veto_dR = 0.1
    # Keep using old TuneP for past selections
    if 'Dec2012' not in Selection.__file__:
        leptons.muon_track_for_momentum = cms.string('TuneP')
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the MuPrescaled paths, we don't care about e-mu events.
        if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
            continue

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)
        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        # Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        if cut_name == 'Simple':
            alldil.electron_cut_mask = cms.uint32(0)
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'#to be changed for first runs
            #alldil.loose_cut = 'isGlobalMuon && pt > 0.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.sort_by_pt = True
            dil.do_remove_overlap = False
            if hasattr(dil, 'back_to_back_cos_angle_min'):
                delattr(dil, 'back_to_back_cos_angle_min')
            if hasattr(dil, 'vertex_chi2_max'):
                delattr(dil, 'vertex_chi2_max')
            if hasattr(dil, 'dpt_over_pt_max'):
                delattr(dil, 'dpt_over_pt_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
        elif 'MuPrescaled' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
            assert alldil.tight_cut == trigger_match
            alldil.tight_cut = prescaled_trigger_match

        # Histos now just needs to know which leptons and dileptons to use.
        histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    pobj = process.muonPhotonMatch * reduce(lambda x,y: x*y, path_list)
    if 'VBTF' not in cut_name and cut_name != 'Simple':
        pobj = process.goodDataFilter * pobj # goodDataFilter is not applied if you run locally the code
#    if 'MuPrescaled' in cut_name: ####### Now it seams that there are no prescaled path ########
#        pobj = process.PrescaleToCommon * pobj ####### Now it seams that there are no prescaled path ########
    path = cms.Path(pobj)
    setattr(process, pathname, path)


def ntuplify(process, fill_gen_info=False):
    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler',
                                           dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlinePrimaryVertices'),
                                           )
    process.SimpleNtuplerEmu = process.SimpleNtupler.clone(dimu_src = cms.InputTag('SimpleMuonsElectronsAllSigns'))

    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        process.SimpleNtupler.hardInteraction = hardInteraction
        
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.SimpleNtupler * process.SimpleNtuplerEmu
ntuplify(process) #to have ntuples also running in interactive way

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.triggerSummaryAnalyzerAOD

    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.pathSimple *= process.PrintOriginalMuons

    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.PrintEventSimple

    #- 2011-2012 selection (Nlayers > 8)
    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu

    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
    if hasattr(process, 'pathOur2012'):
        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
    process.CheckPrescale.trigger_paths = cms.vstring(*trigger_paths)
    process.pCheckPrescale = cms.Path(process.CheckPrescale)
###check_prescale(process, trigger_paths) ### ony to check if it works

def for_data(process):
    process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'
    ntuplify(process)
    #check_prescale(process, trigger_paths) ####### Now it seams that there are no prescaled path ########

def for_mc(process, hlt_process_name, fill_gen_info):
    ntuplify(process, fill_gen_info)
    switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)

def get_dataset(run):
    #JMTBAD common with dataset_details in submit below, make a DataSamples.py?
    run = int(run)
    if 190450 <= run <= 191284:
        return '/SingleMu/tucker-datamc_SingleMuRun2012A_Prompt_190450_191284_20120418134612-57b19813ab8f2ab142c4566dc6738156/USER'
    else:
        raise ValueError('dunno how to do run %i' % run)

if 'int_data' in sys.argv:
    for_data(process)
    printify(process)
    
if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', False)
    printify(process)
    
if 'gogo' in sys.argv:
    for_data(process)
    printify(process)
    
    n = sys.argv.index('gogo')
    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
    print run, lumi, event
    run = int(run)
    lumi = int(lumi)
    event = int(event)
    filename = [x for x in sys.argv if x.endswith('.root')]
    if filename:
        filename = filename[0]
    else:
        dataset = get_dataset(run)
        print dataset
        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
        print repr(output)
        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
    print filename
    process.source.fileNames = [filename]
    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
    set_events_to_process(process, [(run, event)])

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_datamc_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
#config.JobType.priority = 1

config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'phys03'
job_control
config.Data.publication = False
config.Data.publishDataName = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/rradogna'

#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_IT_Legnaro'

'''
    
    just_testing = 'testing' in sys.argv
        
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [
#            ('SingleMuRun2012A_13Jul2012_190450_193751', '/SingleMu/slava-datamc_SingleMuRun2012A-13Jul2012_190450_193751_20121011073628-426a2d966f78bce6bde85f3ed41c07ba/USER'),
#            ('SingleMuRun2012A_06Aug2012_190782_190949', '/SingleMu/slava-datamc_SingleMuRun2012A-recover-06Aug2012_190782_190949_20121011120430-426a2d966f78bce6bde85f3ed41c07ba/USER'),
#            ('SingleMuRun2012B_13Jul2012_193752_196531', '/SingleMu/slava-datamc_SingleMuRun2012B-13Jul2012_193752_196531_20121012044921-426a2d966f78bce6bde85f3ed41c07ba/USER'),

#            ('ExpressPhysicsRun2015B-Express_251161_251252',    '/ExpressPhysics/rradogna-datamc_ExpressPhysicsRun2015B-Express_251161_251252_20150709121231-b1bb0cb9d6986acfdad6df64560b6f76/USER'),
#            ('ExpressPhysicsRun2015B-Express_251253_251521',    '/ExpressPhysics/rradogna-datamc_ExpressPhysicsRun2015B-Express_251253_251521_20150711154626-b1bb0cb9d6986acfdad6df64560b6f76/USER'),
#            ('ExpressPhysicsRun2015B-Express_251522_251562',    '/ExpressPhysics/rradogna-datamc_ExpressPhysicsRun2015B-Express_251522_251562_20150712181215-b1bb0cb9d6986acfdad6df64560b6f76/USER'),
#            ('ExpressPhysicsRun2015B-Express_251563_251643', '/ExpressPhysics/rradogna-datamc_ExpressPhysicsRun2015B-Express_251563_251643_20150714134408-0633e624cb4a2f2023c182912c95c757/USER'),

            ('SingleMuonRun2015B-Prompt_251162_251499',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251162_251499_20150713100409-3aa7688518cb1f1b044caf15b1a9ed05/USER'),
#            ('SingleMuonRun2015B-Prompt_251500_251561',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251500_251561_20150714135629-3aa7688518cb1f1b044caf15b1a9ed05/USER'),
#             ('SingleMuonRun2015B-Prompt_251562_251642',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251562_251642_20150716164649-9996471c14459acaec01707975d1e954/USER'),
#             ('SingleMuonRun2015B-Prompt_251500_251642',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251500_251642_20150717155341-9996471c14459acaec01707975d1e954/USER'),
            ('SingleMuonRun2015B-Prompt_251500_251603',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251500_251603_20150718235715-9996471c14459acaec01707975d1e954/USER'),
            ('SingleMuonRun2015B-Prompt_251613_251883',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251613_251883_20150719000207-9996471c14459acaec01707975d1e954/USER'),
                           
            ('SingleMuonRun2015C-Prompt_253888_254914',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015C-Prompt_253888_254914_20150831150018-681693e882ba0f43234b3b41b1bbc39d/USER'),
                           
            ('SingleMuonRun2015D-Prompt_256629_256842',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015D-Prompt_256629_256842_20150926113604-c9b39dd88dc98b683a1d7cecc8f6c42c/USER'),

            ]

        lumi_lists = [
#            'NoLumiMask'
#            'DCSOnly',
#            'Run2012PlusDCSOnlyMuonsOnly',
            'Run2015MuonsOnly',
#            'Run2015',
            ]

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name

            new_py = open('histos.py').read()
            new_py += "\nfor_data(process)\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 200
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON_MuonPhys_v2.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys.txt'
#config.Data.lumiMask = '%(lumi_mask)s' #######
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
#                os.system('crab submit -c crabConfig.py --dryrun')#debug
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 5000
#config.Data.unitsPerJob  = 100
    ''')

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples

        combine_dy_samples = len([x for x in samples if x.name in []]) > 0
        #combine_dy_samples = len([x for x in samples if x.name in ['dy50','dy100-200', 'dy200-400', 'dy400-500','dy500-700', 'dy700-800']]) > 0
        print 'combine_dy_samples:', combine_dy_samples

        for sample in reversed(samples):
            print sample.name

            new_py = open('histos.py').read()
            sample.fill_gen_info = sample.name in ['dy50','dy100-200', 'dy200-400', 'dy400-500','dy500-700', 'dy700-800']
            new_py += "\nfor_mc(process, hlt_process_name='%(hlt_process_name)s', fill_gen_info=%(fill_gen_info)s)\n" % sample
            
            if sample.is_madgraph:
                print "Madgraph: re-weight applied"
                new_py += "\nHistosFromPAT.useMadgraphWeight= cms.untracked.bool(True)\n"

            if combine_dy_samples and (sample.name == 'zmumu' or 'dy' in sample.name):
                mass_limits = {
#                    'dy50'      : (  50,     100),
                    #'dy120'     : ( 120,     200),
#                    'dy100-200'     : ( 100,     200),
#                    'dy200-400'     : ( 200,     400),
#                    'dy400-500'     : ( 400,     500),
#                    'dy500-700'     : ( 500,     700),
#                    'dy700-800'     : ( 700,     800),
#                    'dy1500'    : (1000,    1500),
#                    'dy2000'    : (1500,    2000),
#                    'dy3000'    : (2000,    3000),
                    #'dy7500'    : (6000,    7500),
                    #'dy8500'    : (8500,    9500),
                    #'dy9500'    : (9500,  100000),
                    }
                lo,hi = mass_limits[sample.name]
                from SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi import dy_gen_mass_cut
                new_cut = dy_gen_mass_cut % locals()

                new_py += '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DYGenMassFilter_cfi')

process.DYGenMassFilter.cut = "%(new_cut)s"
for pn,p in process.paths.items():
    setattr(process, pn, cms.Path(process.DYGenMassFilter*p._seq))
''' % locals()

            open('histos_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % sample)
            if not just_testing:
                os.system('crab submit  -c crabConfig.py') #--dryrun
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

#        if not just_testing:
#            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
