#!/usr/bin/env python

import sys, os, glob, datetime, time
from collections import defaultdict
from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import crab_status, crabify_list, last_crab_dir

just_testing = 'testing' in sys.argv

def do(cmd):
    print cmd
    if not just_testing:
        os.system(cmd)

# Figure out what run to start with, based on the crab dirs that
# already exist.
print '*** finding new runs to use'
last_dir = last_crab_dir()
if 'promptB' in last_dir:
    print 'trying to parse arguments.xml in', last_dir
    runs_done = set()
    for line in open(os.path.join(last_dir, 'share/arguments.xml')):
        # Time for some really crappy, fragile parsing!
        if 'Lumis=' not in line: continue
        line = line.split('Lumis="')[1].split('"')[0]
        for x in line.split(','):
            runs_done.add(int(x.split(':')[0]))
    last_run_done = max(runs_done)
    print 'runs already done in last crab batch: min: %i max: %i' % (min(runs_done), last_run_done)
else:
    last_run_done = 144114 #int(open('last_run_done', 'rt').read())

# Find latest run in dataset.
dataset = '/Mu/Run2010B-PromptReco-v2/RECO'
dbsoutput = os.popen('dbs search --query="find max(run) where dataset=%s"' % dataset).read()
for line in dbsoutput.split('\n'):
    try:
        max_run = int(line)
    except ValueError:
        pass
print 'max run in %s: %i' % (dataset, max_run)

# Find latest here at FNAL.
dbsoutput = os.popen('dbs search --query="find file where dataset=%s"' % dataset).read()
for line in dbsoutput.split('\n'):
    if '.root' in line:
        directory = '/'.join(line.strip().split('/')[:-3])
        run1k = max(os.listdir('/pnfs/cms/WAX/11' + directory))
        directory += '/' + run1k
        run = int(run1k + max(os.listdir('/pnfs/cms/WAX/11' + directory)))
        if run < max_run:
            print 'max run here at fnal: %i' % run
            max_run = run
        break

print 'max_run: %i  last done: %i' % (max_run, last_run_done)
if max_run <= last_run_done:
    print 'no new runs, quitting.'
    sys.exit(0)

print '*** submitting jobs'
do('python tuple_data.py submit %i %i' % (last_run_done + 1, max_run))

# Spin until the crab jobs are all done.
last_dir = last_crab_dir()
print '*** spinning on', last_dir
while not just_testing:
    status = crab_status(last_dir)
    if 'Done' in status.keys():
        os.system('crab -c %s -getoutput %s' % (last_dir, crabify_list(status['Done'])))
        status = crab_status(last_dir)
        
    print 'jobs status at', datetime.datetime.now().ctime(), status

    if 'Submitted' in status.keys() or 'Ready' in status.keys() or 'Running' in status.keys():
        print 'sleeping fifteen minutes...'
        time.sleep(15*60)
        continue

    if status.keys() != ['Retrieved_0_0']:
        raise RuntimeError('some problem with the crab jobs for %s !' % last_dir)

    break

do('crab -c %s -report' % last_dir)

# Merge the output into one PAT tuple file.
print '*** crab jobs done! merging output:'
do('cmsRun merge.py %s' % last_dir)

print '*** merging done. Making histos.'
in_fn = os.path.join(last_dir, 'res', 'merged.root')
out_fn_base = os.path.join('ana_datamc', 'ana_' + last_dir.replace('crab/crab_', '').replace('datamc', 'datamc_data') + '_%s.rout')
for x in ['allgood', 'zp2mu']:
    lumis_fn = '/afs/cern.ch/user/t/tucker/runreg/output/full_%s.cmssw' % x
    json_fn = lumis_fn.replace('.cmssw', '.json')
    out_fn = out_fn_base % x
    for_lumi_fn = out_fn.replace('.rout', '.forlumi.json')
    do('cp %s %s' % (lumis_fn, out_fn.replace('.rout', '.cmssw')))
    do('cp %s %s' % (json_fn,  out_fn.replace('.rout', '.runreg.json')))
    do('compareJSON.py --and %s %s > %s' % (os.path.join(last_dir, 'res', 'lumiSummary.json'), json_fn, for_lumi_fn))
    do('lumiCalc.py -i %s overview > %s' % (for_lumi_fn, out_fn.replace('.rout', '.lumi')))
    do('cmsRun histos.py data %s %s %s' % (in_fn, lumis_fn, out_fn))

print '** done!'
