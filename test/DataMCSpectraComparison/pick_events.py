#!/usr/bin/env python

# The official edmPickEvents tool is unwieldy, especially having to
# launch across two or more datasets using crab, etc. So, here's this.

import sys, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process

process = cms.Process('PickEvent')
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('dummy.root'))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.out = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pickevents.root'))
process.end = cms.EndPath(process.out)

if __name__ == '__main__' and 'submit' in sys.argv:
    import datetime, os
    from collections import defaultdict
    from FWCore.PythonUtilities.LumiList import LumiList
    from pprint import pformat
    
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '%(name)s' 
config.General.workArea = 'Pickevents_%(uniq)s_%(name)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'pick_events_crab.py'   
config.JobType.priority = 1

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 10000
config.Data.publication = True
config.Data.publishDataName = '%(name)s'
config.Data.outLFN = '/store/user/federica/Pickevents' 

config.Site.storageSite = 'T2_US_Purdue'

'''

    just_testing = 'testing' in sys.argv
    create_only = 'create_only' in sys.argv
    uniq = datetime.datetime.today().strftime('%Y%m%d%H%M%S')

    def get_dataset(run):
        if 165071 <= run <= 168437:
            return '/SingleMu/Run2011A-PromptReco-v4/RECO'
        elif 170053 <= run <= 172619:
            return '/SingleMu/Run2011A-PromptReco-v5/RECO'
        elif 172620 <= run <= 175770:
            return '/SingleMu/Run2011A-PromptReco-v6/RECO'
        elif 175832 <= run <= 180296:
            return '/SingleMu/Run2011B-PromptReco-v1/RECO'
        elif 160329 <= run <= 163869:
            return '/SingleMu/Run2011A-May10ReReco-v1/RECO'
        else:
            raise ValueError('dunno how to do run %i' % run)

    def name_for_dataset(dataset):
        return dataset.replace('/', '').replace('-', '')
    
    batches = defaultdict(list)
    
    # The format of the input file must be one or more lines of
    #   run:lumi:event[:dataset]
    # where 'dataset' is an optional override for that particular
    # event for the dataset to use. Otherwise, we figure out the
    # dataset here based on the run number.
    events_fn = [x for x in sys.argv[1:] if os.path.isfile(x)][0]
    is_mc = False
    for line in open(events_fn):
        line = line.strip()
        if not line:
            continue
        line = line.split(':')
        dataset = None
        if len(line) == 4:
            line, dataset = line[:3], line[3]
        assert len(line) == 3
        run, lumi, event = line
        run, lumi, event = int(run), int(lumi), int(event)
        if dataset is None:
            dataset = get_dataset(run)
        batches[dataset].append((run, lumi, event))
    
    for dataset, batch in batches.iteritems():
        name = name_for_dataset(dataset)
        print name
            
        lumi_mask = []
        events_to_process = []

        for run, lumi, event in batch:
            lumi_mask.append((run, lumi))
            events_to_process.append((run, event))

        ls = set(l for r,l in lumi_mask)

        if ls == set([-1]):
            is_mc = True
        elif -1 in ls:
            raise ValueError('batch for dataset %s has lumis -1 and others' % dataset)
        else:
            is_mc = False

        if not is_mc:
            job_control = '''
lumi_mask = pick_events.json
total_number_of_lumis = -1
lumis_per_job = 1'''
            ll = LumiList(lumis=lumi_mask)
            ll.writeJSON('pick_events.json')
        else:
            job_control = '''
total_number_of_events = -1
events_per_job = 100000'''

        scheduler = 'condor' if 'condor' in sys.argv else 'glite'
        open('crab.py', 'wt').write(crab_cfg % locals())

        pset = open('pick_events.py').read()
        pset += '\nevents_to_process = '
        pset += pformat(events_to_process)
        pset += '\nset_events_to_process(process, events_to_process)\n'
        open('pick_events_crab.py', 'wt').write(pset)

        if not just_testing:
            if create_only:
                os.system('crab submit')
        else:
            print 'crab.py'
            os.system('less crab.py')
            print 'pick_events.json'
            os.system('less pick_events.json')
            print 'pick_events_crab.py'
            os.system('diff pick_events.py pick_events_crab.py | less')

    if not just_testing:
        os.system('rm crab.py pick_events.json pick_events_crab.py pick_events_crab.pyc')
