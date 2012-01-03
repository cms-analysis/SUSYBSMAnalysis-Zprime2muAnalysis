#!/usr/bin/env python

import sys, os, datetime, FWCore.ParameterSet.Config as cms
from tuple_common import process, crab_cfg

process.source.fileNames = ['/store/data/Run2011A/SingleMu/AOD/08Nov2011-v1/0003/FC923030-7413-E111-9584-1CC1DE1CEDB2.root']
process.GlobalTag.globaltag = 'FT_R_44_V11::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and 'submit' in sys.argv:
    job_control_ex = '''
total_number_of_lumis = -1
lumis_per_job = %(lumis_per_job)s
%(lumi_mask)s
'''

    lumis_per_job = 300
    lumi_mask = ''

    create_only = 'create_only' in sys.argv
    just_testing = 'testing' in sys.argv
    scheduler = 'condor' if 'grid' not in sys.argv else 'glite'

    def submit(d):
        new_py = open('tuple_data.py').read()
        new_py += '\n\nprocess.GlobalTag.globaltag = "%(tag)s::All"\n' % d
        pset = 'crab/psets/tuple_data_crab_%(name)s.py' % d
        open(pset, 'wt').write(new_py)

        job_control = job_control_ex % d
        for k,v in locals().iteritems():
            d[k] = v
        open('crab.cfg', 'wt').write(crab_cfg % d)
        if not just_testing:
            if create_only:
                os.system('crab -create')
            else:
                os.system('crab -create -submit all')
            os.system('rm -f crab.cfg tmp.json')

    run_limits = []
    for x in sys.argv:
        try:
            run_limits.append(int(x))
        except ValueError:
            pass

    if run_limits:
        raise NotImplementedError('no run_limits in 44X')
    
        run1, run2 = run_limits
        if len(run_limits) != 2 or run1 > run2:
            raise RuntimeError('if any, must specify exactly two numeric arguments   min_run max_run  with max_run >= min_run')

        # Make up a fake lumi_mask that contains all lumis possible
        # for every run in the run range, since crab doesn't seem to
        # listen for a runselection parameter anymore.
        json = ['"%i": [[1,26296]]' % r for r in xrange(run_limits[0], run_limits[1] + 1)]
        open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')
        lumi_mask = 'lumi_mask = tmp.json'

        name = 'SingleMuRun2011B_Prompt_%i_%i_%s' % (run_limits[0], run_limits[1], datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        print name

        if run1 >= 999999:
            dataset = '/Nope/NoDataset/NO'
        else:
            raise ValueError("don't know how to do a run_limits production for run range [%i,%i]" % run_limits)
        
        tag = 'FT_R_44_V9'
        submit(locals())
    else:
        x = [
            ('SingleMuRun2011A_Nov08', '/SingleMu/Run2011A-08Nov2011-v1/AOD', 'FT_R_44_V9'),
            ('SingleMuRun2011B_Nov19', '/SingleMu/Run2011B-19Nov2011-v1/AOD', 'FT_R_44_V11'),
            ]
        for name, dataset, tag in x:
            submit(locals())

