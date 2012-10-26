#!/usr/bin/env python

# This script dumps all the L1 and HLT prescales from the DB for a
# given HLT path (e.g. HLT_Mu30_v5) and alerts if any prescale is not
# 1.

import sys, os, re
from collections import defaultdict
from RecoLuminosity.LumiDB import lumiQueryAPI
from FWCore.PythonUtilities.LumiList import LumiList

# Bomb if args not specified.
if len(sys.argv) < 3 or not os.path.isfile(sys.argv[1]):
    print 'usage: prescales.py lumis.json pathsubstr'
    sys.exit(1)
json = sys.argv[1]
path = sys.argv[2]

ll = LumiList(json)

# Why doesn't LumiList have a getRunsAndLumis method? It can be
# constructed from a dict of runsAndLumis... Anyway we'll use sets for
# faster searching later.
runs_and_lumis = defaultdict(set)
for r,l in ll.getLumis():
    runs_and_lumis[r].add(l)
runs = runs_and_lumis.keys()
runs.sort()

# The path to look for: e.g. HLT_Mu30_v*. Should only find one per
# run, i.e. v4 and v5 do not exist simultaneously.
path_re = re.compile(r'HLT_%s_v\d+' % path)

# Magic.
parameters = lumiQueryAPI.ParametersObject()
session, svc = lumiQueryAPI.setupSession('frontier://LumiCalc/CMS_LUMI_PROD', None, parameters, False)
session.transaction().start(True)
schema = session.nominalSchema()

# Loop over all the requested runs, using the selected lumis, and get
# the prescales.
l1_prescales = {}
hlt_prescales = {}
for run in runs:
    lumis = runs_and_lumis[run]

    # Find the trigger path requested, including its L1 seed. If we
    # don't find exactly one satisfying the RE, then raise an
    # error. If the L1 seed is the AND or OR of any triggers, also
    # raise an error.
    q = schema.newQuery()
    triggers = lumiQueryAPI.hlttrgMappingByrun(q, run)
    del q
    if not triggers:
        print 'Run %i: warning, no trigger mapping available!' % run
        continue
    paths = []
    for hlt_path, l1_seed in triggers.iteritems():
        mo = path_re.search(hlt_path)
        if mo is not None:
            paths.append((hlt_path, l1_seed))
    if len(paths) != 1:
        raise RuntimeError('for run %i, did not find exactly one path matching the specification: paths = %s' % (run, repr(paths)))
    hlt_path, l1_seed = paths[0]
    l1_seed_l = l1_seed.replace('"', '').replace(' AND ', 'SPLITHERE').replace(' OR ', 'SPLITHERE').split('SPLITHERE')
    if len(l1_seed_l) != 1:
        raise RuntimeError('for run %i, not exactly one L1 seed for HLT path %s: %s %s' % (run, hlt_path, l1_seed, repr(l1_seed_l)))
    l1_seed = l1_seed_l[0]

    # Get the HLT prescales for the above path.
    q = schema.newQuery()
    hlt_info = lumiQueryAPI.hltBypathByrun(q, run, hlt_path)
    del q
    if not hlt_info:
        print 'Run %i: warning, no hlt_info available for path %s!' % (run, hlt_path)
        continue
    hlt_prescales_seen = set()
    for ls, (input, accept, prescale) in hlt_info.iteritems():
        if ls in lumis:
            hlt_prescales_seen.add(prescale)
            hlt_prescales[(run,ls)] = prescale

    # Get the L1 prescale too.
    q = schema.newQuery()
    l1_info = lumiQueryAPI.trgBybitnameByrun(q, run, l1_seed)
    del q
    l1_prescales_seen = set()
    for ls, (count, deadtime, bit_number, prescale) in l1_info.iteritems():
        if ls in lumis:
            l1_prescales_seen.add(prescale)
            l1_prescales[(run,ls)] = prescale

    print 'Run %i:  HLT path, prescales: %s %s   L1 seed, prescales: %s %s' % (run, hlt_path, repr(sorted(hlt_prescales_seen)), l1_seed, repr(sorted(l1_prescales_seen)))

session.transaction().commit()
del session
del svc

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import to_pickle
to_pickle((l1_prescales, hlt_prescales), path + '.gzpickle')

'''
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import from_pickle
l1, hlt = from_pickle('Mu15.gzpickle')
rls = sorted(set(l1.keys()) & set(hlt.keys()))
max_seen = 0
for rl in rls:
    prescale = l1[rl]*hlt[rl]
    if prescale > max_seen:
        print 'new max prescale seen for', rl, ':', prescale
        max_seen = prescale
'''
