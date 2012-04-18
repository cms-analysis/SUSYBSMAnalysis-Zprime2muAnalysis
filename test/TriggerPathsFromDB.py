#!/usr/bin/env python

import os, re, sys
from collections import defaultdict
from itertools import combinations
from FWCore.PythonUtilities.LumiList import LumiList

# Should rewrite not to hit the db for every cfg, but just get the HLT
# key for each run and then only get cfgs for unique keys.

dcsonly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/DCSOnly/json_DCSONLY.txt') # JMTBAD use DCSOnly_ll from goodlumis once sorted
runs = sorted(int(run) for run in dcsonly_ll.getRuns())
cmd = 'edmConfigFromDB --cff --runNumber %i --noedsources --noes --noservices --nomodules'

path_re           = re.compile(r'(HLT_Mu40_eta2p1_v\d+)')
prescaled_path_re = re.compile(r'(HLT_Mu15_eta2p1_v\d+)')

paths_and_filters = defaultdict(list)

for run in runs:
    print 'run:', run,
    sys.stdout.flush()
    path = prescaled_path = filter = prescaled_filter = None
    for line in os.popen(cmd % run):
        if 'cms.Path' not in line:
            continue
        filt = line.split(' + ')[-2] # JMTBAD fragile
        mo = path_re.search(line)
        if mo is not None:
            path = mo.group(1)
            filter = filt
        mo = prescaled_path_re.search(line)
        if mo is not None:
            prescaled_path = mo.group(1)
            prescaled_filter = filt
    pandf = (path, filter, prescaled_path, prescaled_filter)
    if None in pandf and set(pandf) != set([None]):
        raise ValueError('for run %i, one of the paths/filters was not found: path: %s  filter: %s  prescaled_path: %s  prescaled_filter: %s (either all must be found or none may be found)' % (run, path, filter, prescaled_path, prescaled_filter))
    if None not in pandf:
        paths_and_filters[pandf].append(run)
    print '%20s %70s %20s %70s' % pandf

print

# Before printing out simple summary lines for copy/paste into
# Zprime2muTriggerPathsAndFilters, check to make sure the simple way
# makes sense. This means no overlapping run ranges are allowed.

run_ranges = paths_and_filters.values()
for rr1, rr2 in combinations(run_ranges, 2):
    rr1 = set(xrange(min(rr1), max(rr1)+1))
    rr2 = set(xrange(min(rr2), max(rr2)+1))
    if rr1 & rr2:
        raise ValueError('overlap') # JMTBAD fix error message

print '\nno overlap found, so can keep simple. to paste (resorting up to you):'
code_line = 'else if (run >= %(run1)i && run <= %(run2)i) { path = "%(path)s", filter = "%(filter)s", prescaled_path = "%(prescaled_path)s", prescaled_filter = "%(prescaled_filter)s"; }'
for pandf in sorted(paths_and_filters.iterkeys()):
    runs = paths_and_filters[pandf]
    run1, run2 = min(runs), max(runs)
    path, filter, prescaled_path, prescaled_filter = pandf
    print code_line % locals()

