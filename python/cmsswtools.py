#!/usr/bin/env python

import sys, FWCore.ParameterSet.Config as cms

def files_from_argv(process):
    files = [x for x in sys.argv if '.root' in x]
    for f in sys.argv:
        if '.root' not in x:
            continue
        
        if '/store' in f:
            files[i] = f[f.index('/store'):]
        else:
            files[i] = 'file:' + f
        
    process.source.fileNames = files
    
def set_events_to_process(process, run_events, run=None):
    '''Set the PoolSource parameter eventsToProcess appropriately,
    given the desired runs/event numbers passed in. If run is None,
    run_events must be a list of 2-tuples, each entry being a (run,
    event) pair. Otherwise, the run number is taken from run, and
    run_events is just a list of event numbers to be paired with run.

    run_events can also be a list of 3-tuples, where the middle entry
    in each is the lumisection number. This is ignored for now.
    '''
    if run is not None:
        run_events = [(run, event) for event in run_events]
    #er = cms.untracked.VEventID(*[cms.untracked.EventID(*x) for x in run_events])
    er = cms.untracked.VEventRange(*[cms.untracked.EventRange(x[0],x[-1],x[0],x[-1]) for x in run_events])
    process.source.eventsToProcess = er        
