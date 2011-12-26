#!/usr/bin/env python

import os, sys
from collections import defaultdict
from pprint import pprint
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, ttree_iterator

f = ROOT.TFile(sys.argv[1])
t = f.Mu15.Get('t')
print t.GetEntriesFast(), 'entries in tree'

de = defaultdict(list)
dp = defaultdict(list)
for jentry, tt in ttree_iterator(t):
    if jentry % 2000000 == 0:
        print jentry
    k = (t.run, t.lumi)
    de[k].append(t.event)
    dp[k].append((t.l1,t.hlt))

for k,v in de.iteritems():
    assert len(v) == len(set(v))

for k,v in dp.iteritems():
    if len(set(v)) != 1:
        print 'run/LS %s with %i events has multiple prescale values: %s' % (k, len(de[k]), sorted(set(v)))

lscounts = defaultdict(int)
evcounts = defaultdict(int)
uniq = set()
for v in dp.itervalues():
    uniq.update(v)
    for l1hlt in set(v):
        lscounts[l1hlt] += 1
    for l1hlt in v:
        evcounts[l1hlt] += 1

uniq = sorted(set(uniq))
print 'unique prescale factors'
pprint(uniq)
print 'event counts'
pprint(dict(evcounts))
print 'ls counts'
pprint(dict(lscounts))



