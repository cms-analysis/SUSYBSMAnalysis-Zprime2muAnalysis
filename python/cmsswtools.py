#!/usr/bin/env python

import glob, os, sys, FWCore.ParameterSet.Config as cms

def cmssw_version(as_tuple=True):
    if not os.environ.has_key('CMSSW_VERSION'):
        raise RuntimeError('CMSSW_VERSION environment variable not set!')
    x = os.environ['CMSSW_VERSION'].replace('CMSSW_', '')
    if as_tuple:
        x = x.split('_')
        return tuple(int(y) for y in x[:3]) + tuple(x[3:])
    else:
        return x
                                
def files_from_argv(process):
    files = []
    for f in sys.argv:
        if '.root' not in f:
            continue
        if '/store' in f:
            f = f[f.index('/store'):]
        else:
            f = 'file:' + f
        files.append(f)
    if not files:
        raise RuntimeError('no .root files found in argv: %s' % repr(sys.argv))
    process.source.fileNames = files

def files_from_path(process, path):
    files = ['file:%s' % x for x in glob.glob(os.path.join(path, '*.root'))]
    if not files:
        raise RuntimeError('no .root files found in %s' % path)
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

def set_preferred_alignment(process, name, connect, **kwargs):
    '''Function to select a set of alignment constants. Useful when doing
    track re-reconstruction.

    Example use:

    extra_alignment = [('frontier://FrontierProd/CMS_COND_31X_FROM21X', {'CSCAlignmentRcd': 'CSCAlignmentRcd_CRAFT_PG-hardware-globalMuons_v3_offline'})]
    for i, (connect, rcds) in enumerate(extra_alignment):
        set_preferred_alignment(process, 'extraAlignment%i' % i, connect, **rcds)
    '''

    if len(kwargs) == 0:
        raise ValueError, 'must specify at least one record parameter'
    
    from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
    alignment_source = cms.ESSource('PoolDBESSource',
        CondDBSetup,
        connect = cms.string(connect),
        toGet = cms.VPSet()
    )

    for record, tag in kwargs.iteritems():
        alignment_source.toGet.append(cms.PSet(record = cms.string(record), tag = cms.string(tag)))

    setattr(process, name, alignment_source)
    process.prefer(name)

    #code = "cms.ESPrefer('PoolDBESSource', name, %s)" % ', '.join(['%s=cms.vstring("%s")' % x for x in kwargs.iteritems()])
    #setattr(process, name + '_es_prefer', eval(code))
