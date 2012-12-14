#!/usr/bin/env python

import sys, os, datetime, FWCore.ParameterSet.Config as cms
from tuple_common import process, crab_cfg

process.source.fileNames = ['/store/data/Run2012A/SingleMu/AOD/13Jul2012-v1/00000/009C369E-85D0-E111-BD58-1CC1DE046FC0.root']
process.GlobalTag.globaltag = 'FT_53_V10_AN2::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    job_control_ex = '''
total_number_of_lumis = -1
lumis_per_job = %(lumis_per_job)s
%(lumi_mask)s
'''

    lumis_per_job = 200
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
        run1, run2 = run_limits
        if len(run_limits) != 2 or run1 > run2:
            raise RuntimeError('if any, must specify exactly two numeric arguments   min_run max_run  with max_run >= min_run')

        # Make up a fake lumi_mask that contains all lumis possible
        # for every run in the run range, since crab doesn't seem to
        # listen for a runselection parameter anymore.
        json = ['"%i": [[1,26296]]' % r for r in xrange(run_limits[0], run_limits[1] + 1)]
        open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')
        lumi_mask = 'lumi_mask = tmp.json'

        if run1 == 190782 and run2 == 190949:
            # Special settings for 6-Aug reprocessing of 5 runs
            dataset = '/SingleMu/Run2012A-recover-06Aug2012-v1/AOD'
            name    = 'SingleMuRun2012A-recover-06Aug2012'
            tag     = 'FT_53_V10_AN2'
        elif run1 >= 190450 and run1 < 193752:
            dataset = '/SingleMu/Run2012A-13Jul2012-v1/AOD'
            name    = 'SingleMuRun2012A-13Jul2012'
            tag     = 'FT_53_V10_AN2'
        elif run1 >= 193752 and run1 < 196532:
            dataset = '/SingleMu/Run2012B-13Jul2012-v1/AOD'
            name    = 'SingleMuRun2012B-13Jul2012'
            tag     = 'FT_53_V10_AN2'
        elif run1 >= 197556 and run1 < 198914:
            dataset = '/SingleMu/Run2012C-24Aug2012-v1/AOD'
            name    = 'SingleMuRun2012C-24Aug2012'
            tag     = 'FT_53_V10_AN2'
        elif run1 >= 198934 and run1 < 203773:
            dataset = '/SingleMu/Run2012C-PromptReco-v2/AOD'
            name    = 'SingleMuRun2012C-Prompt'
            tag     = 'GR_P_V42_AN2'
        elif run1 >= 203773:
            dataset = '/SingleMu/Run2012D-PromptReco-v1/AOD'
            name    = 'SingleMuRun2012D-Prompt'
            tag     = 'GR_P_V42_AN2'
        else:
            raise ValueError("don't know how to do a run_limits production for run range [%i,%i]" % run_limits)

        name = '%s_%i_%i_%s' % (name, run_limits[0], run_limits[1], datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        print name, tag

        submit(locals())
    else:
        raise ValueError('must do a run-limits production until one dataset is closed')
        x = [
            ]
        for name, dataset, tag in x:
            submit(locals())

