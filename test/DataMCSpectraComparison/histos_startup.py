#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

process.source.fileNames =[#'file:./pat.root',

#                           '/store/user/rradogna/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/datamc_dy50/150708_150131/0000/pat_110.root',

                           '/store/user/rradogna/SingleMuon/datamc_SingleMuonRun2015C-Prompt_253888_254914_20150831150018/150831_130042/0000/pat_452.root',

                           ]
process.maxEvents.input = -1
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'## solo per proare i dati
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A'## solo per proare i mc
#process.GlobalTag.globaltag = '74X_mcRun2_startup_realistic50ns_v0' #mc startup
#process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1 # default 1000

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
#HistosFromPAT.useMadgraphWeight= True #if we want to re-wigth madgraph samples. For other samples the weight is always 1. Does not work with data.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi')
process.HistosFromPAT.leptonsFromDileptons = True
# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff as OurSelectionDec2012


from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold

##################
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels, rec_level_module
tracks = ['inner', 'tunepnew','startup']
rec_levels(process, tracks)
##################
#Zprime2muAnalysis_cfg imports #Zprime2muAnalysis_cfg where is defined Zprime2muAnalysisSequence = cms.Sequence(muonPhotonMatch * leptons * allDimuons * dimuons)
# full sequence process.p = cms.Path(process.goodDataFilter*process.muonPhotonMatch * process.leptons * process.allDimuons * process.dimuons * process.HistosFromPAT)



process.p = cms.Path(process.goodDataFilter*process.Zprime2muAnalysisSequence)

process.p *= rec_level_module(process, process.HistosFromPAT, 'Histos', tracks)


def rec_level_module_ntupler(process, module, name, tracks):
    p = []
    for t in tracks:
        h = module.clone()
        if hasattr(h, 'dimu_src'):
            h.dimu_src = cms.InputTag('dimuons' + t)
        setattr(process, name + t, h)
        p.append(h)
    return reduce(lambda x,y: x*y, p)

def ntuplify(process, fill_gen_info=False):
    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler',
                                           dimu_src = cms.InputTag('dimuons'),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlinePrimaryVertices'),
                                           )
    print fill_gen_info
    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        process.SimpleNtupler.hardInteraction = hardInteraction
    process.p *= rec_level_module_ntupler(process, process.SimpleNtupler, 'SimpleNtupler', tracks)

#ntuplify(process, True) #to have ntuples also running in interactive way

#def printify(process):
#    process.MessageLogger.categories.append('PrintEvent')
#
#    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
#    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
#    if hasattr(process, 'pathSimple'):
#        process.pathSimple *= process.triggerSummaryAnalyzerAOD
#
#    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
#    process.pathSimple *= process.PrintOriginalMuons
#
#    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
#    if hasattr(process, 'pathSimple'):
#        process.pathSimple *= process.PrintEventSimple

#    #- 2011-2012 selection (Nlayers > 8)
#    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
#    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
#    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
#    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu
#
#    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
#    if hasattr(process, 'pathOur2012'):
#        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
#        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
#        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
#        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu


def for_data(process):
    process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'
    ntuplify(process)


def for_mc(process, hlt_process_name, fill_gen_info):
    ntuplify(process, fill_gen_info)



if 'int_data' in sys.argv:
    for_data(process)
#    printify(process)

if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', False)
#    printify(process)

#if 'gogo' in sys.argv:
#    for_data(process)
##    printify(process)
#
#    n = sys.argv.index('gogo')
#    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
#    print run, lumi, event
#    run = int(run)
#    lumi = int(lumi)
#    event = int(event)
#    filename = [x for x in sys.argv if x.endswith('.root')]
#    if filename:
#        filename = filename[0]
#    else:
#        dataset = get_dataset(run)
#        print dataset
#        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
#        print repr(output)
#        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
#    print filename
#    process.source.fileNames = [filename]
#    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
#    set_events_to_process(process, [(run, event)])



if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_datamc_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_startup_crab.py'
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


            ('SingleMuonRun2015B-Prompt_251162_251499',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251162_251499_20150713100409-3aa7688518cb1f1b044caf15b1a9ed05/USER'),

            ('SingleMuonRun2015B-Prompt_251500_251603',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251500_251603_20150718235715-9996471c14459acaec01707975d1e954/USER'),
            ('SingleMuonRun2015B-Prompt_251613_251883',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015B-Prompt_251613_251883_20150719000207-9996471c14459acaec01707975d1e954/USER'),
            
            ('SingleMuonRun2015C-Prompt_253888_254914',    '/SingleMuon/rradogna-datamc_SingleMuonRun2015C-Prompt_253888_254914_20150831150018-681693e882ba0f43234b3b41b1bbc39d/USER'),

            ]

        lumi_lists = [
#            'NoLumiMask'
#            'DCSOnly',
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

            new_py = open('histos_startup.py').read()
            new_py += "\nfor_data(process)\n"
            open('histos_startup_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 50
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.lumiMask = '%(lumi_mask)s' #######
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py --dryrun')
#                os.system('crab submit -c crabConfig.py --dryrun')#debug
            else:
                cmd = 'diff histos_startup.py histos_startup_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
            os.system('rm crabConfig.py histos_startup_crab.py histos_startup_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 10000
#config.Data.unitsPerJob  = 100
    ''')

        from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples

        combine_dy_samples = len([x for x in samples if x.name in []]) > 0
#        combine_dy_samples = len([x for x in samples if x.name in ['dy50','dy100to200', 'dy200to400', 'dy400to500','dy500to700', 'dy700to800', 'dy800to1000']]) > 0
        print 'combine_dy_samples:', combine_dy_samples

        for sample in reversed(samples):
            print sample.name

            new_py = open('histos_startup.py').read()
            sample.fill_gen_info = sample.name in ['dy50','dy100to200']
            new_py += "\nfor_mc(process, hlt_process_name='%(hlt_process_name)s', fill_gen_info=%(fill_gen_info)s)\n" % sample
            #### NEW: manage madgraph weights
            if sample.is_madgraph:
                print "Madgraph sample: weights will be applied"
                new_py += "\nHistosFromPAT.useMadgraphWeight= True\n"
                new_py += "\nprocess.Histosinner.useMadgraphWeight= True\n"
                new_py += "\nprocess.Histostunepnew.useMadgraphWeight= True\n"
                new_py += "\nprocess.Histosstartup.useMadgraphWeight= True\n"
            
            if combine_dy_samples and (sample.name == 'zmumu' or 'dy' in sample.name):
                mass_limits = {
                    'dy50'      : (  50,     100),
                    #'dy120'     : ( 120,     200),
#                    'dy100to200'     : ( 100,     200),
#                    'dy200to400'     : ( 200,     400),
#                    'dy400to500'     : ( 400,     500),
#                    'dy500to700'     : ( 500,     700),
#                    'dy700to800'     : ( 700,     800),
#                    'dy800to1000'     : ( 800,     1000),
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

            open('histos_startup_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % sample)
            if not just_testing:
#                os.system('crab submit --dryrun -c crabConfig.py')
                os.system('crab submit -c crabConfig.py --dryrun')
            else:
                cmd = 'diff histos.py histos_startup_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            os.system('rm crabConfig.py histos_startup_crab.py histos_startup_crab.pyc')

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()