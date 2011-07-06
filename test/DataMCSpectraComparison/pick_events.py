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
[CMSSW]
%(job_control)s
pset = pick_events_crab.py
datasetpath = %(dataset)s
get_edm_output = 1

[USER]
return_data = 1
ui_working_dir = crab/crab_pickevents_%(uniq)s_%(name)s

[CRAB]
scheduler = %(scheduler)s
jobtype = cmssw
'''

    just_testing = 'testing' in sys.argv
    create_only = 'create_only' in sys.argv
    uniq = datetime.datetime.today().strftime('%Y%m%d%H%M%S')

    def get_dataset(run):
        if run >= 165071:
            return '/SingleMu/Run2011A-PromptReco-v4/AOD'
        return '/SingleMu/Run2011A-May10ReReco-v1/AOD'

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
        open('crab.cfg', 'wt').write(crab_cfg % locals())

        pset = open('pick_events.py').read()
        pset += '\nevents_to_process = '
        pset += pformat(events_to_process)
        pset += '\nset_events_to_process(process, events_to_process)\n'
        open('pick_events_crab.py', 'wt').write(pset)

        if not just_testing:
            if create_only:
                os.system('crab -create')
            else:
                os.system('crab -create -submit')
        else:
            print 'crab.cfg'
            os.system('less crab.cfg')
            print 'pick_events.json'
            os.system('less pick_events.json')
            print 'pick_events_crab.py'
            os.system('diff pick_events.py pick_events_crab.py | less')

    if not just_testing:
        os.system('rm crab.cfg pick_events.json pick_events_crab.py pick_events_crab.pyc')
