#!/usr/bin/env python

from collections import defaultdict
from pprint import pprint
from FWCore.PythonUtilities.LumiList import LumiList
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, ttree_iterator

me = defaultdict(list)
vbtf = defaultdict(list)

me_json = LumiList('sept17.forlumi.json')
vbtf_json = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON.txt')

peak_only = False

#fns = ['ana_datamc_vbtfvsours/ana_datamc_data_jul15_prompt.root', 'ana_datamc_vbtfvsours/ana_datamc_data_promptB_allgood.root']
#fns = ['ana_datamc_data.root', 'ana_datamc_vbtfvsours/ana_datamc_data_promptB_allgood.root']
fns = ['sept17.root']
for fn in fns:
    f = ROOT.TFile(fn)
    for j,t in ttree_iterator(f.SimpleNtuplerVBTF.Get('t')):
        if not peak_only or (t.dil_mass > 60 and t.dil_mass < 120):
            r,l,e,m = t.run, t.lumi, t.event, t.dil_mass
            for x in 'rle':
                exec '%s = int(%s)' % (x,x) # python ints are c longs
            me[(r,l,e)].append(m)

for line in open('log_goldenZmumuCands_NoMassCuts_34p86.log'):
    if '*' not in line or 'Row' in line:
        continue

    line = line.split('*')
    if len(line) != 8:
        continue
    
    r,l,e,m = [x.strip() for x in line[3:7]]
    r,l,m = int(r), int(l), float(m)
    
    # sometimes the number of digits in the event number is too big
    # for the column and root decides to put it in scientific
    # notation...
    try:
        e = int(e)
    except ValueError:
        print 'warning, had to get event number from sci. not.'
        if me.has_key((r,l)):
            absmax = 0.01
            #if len(me[(r,l)]) == 1:
            #    absmax = 1
            for me_e, me_m in me[(r,l)]:
                if abs(me_m - m) < absmax:
                    e = me_e
            if type(e) != long:
                print 'couldnt match', r,l,e,m,me[(r,l)]

    vbtf[(r,l,e)].append(m)

mek = set(me.keys())
vbtfk = set(vbtf.keys())

def m2s(l):
    return ' * '.join('%.1f' % m for m in l)

print 'counts:\tme\tvbtf'
print 'events:\t%i\t%i' % (len(mek), len(vbtfk))
print
print 'diffs in dimu counts?'
print '%6s%6s%14s%14s%14s%25s%25s' % ('mjs','vjs','run','lumi','event', 'me', 'vbtf')
for rle in sorted(mek & vbtfk):
    if len(me[rle]) != len(vbtf[rle]):
        r,l,e = rle
        print '%6s%6s%14s%14s%14s%25s%25s' % (me_json.contains((r,l)), vbtf_json.contains((r,l)), r,l,e, m2s(me[rle]), m2s(vbtf[rle]))
print
print 'symmetric difference:'
print '%6s%6s%14s%14s%14s%25s%25s' % ('mjs','vjs','run','lumi','event', 'me', 'vbtf')
for rle in sorted(mek.symmetric_difference(vbtfk)):
    r,l,e = rle
    print '%6s%6s%14s%14s%14s%25s%25s' % (me_json.contains((r,l)), vbtf_json.contains((r,l)), r,l,e, m2s(me[rle]), m2s(vbtf[rle]))
