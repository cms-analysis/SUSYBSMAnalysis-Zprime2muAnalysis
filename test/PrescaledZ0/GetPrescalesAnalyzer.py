#!/usr/bin/env python

import os, sys
from collections import defaultdict
from pprint import pprint
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, detree

print sys.argv

f = ROOT.TFile(sys.argv[1])
t = f.Get('%s/t' % sys.argv[2])
print t.GetEntriesFast(), 'entries in tree'

de = defaultdict(list)
dp = defaultdict(list)

for (run, lumi, event, l1, hlt) in detree(t, 'run:lumi:event:l1:hlt'):
    de[(run, lumi)].append(event)
    dp[(run, lumi)].append((l1,hlt))

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



