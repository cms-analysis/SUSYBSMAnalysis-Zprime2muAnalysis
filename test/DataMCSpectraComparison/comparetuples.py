#!/usr/bin/env python

from collections import defaultdict
from FWCore.PythonUtilities.LumiList import LumiList
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, ttree_iterator

def doit(fn, tn='SimpleNtupler'):
    me = defaultdict(list)
    f = ROOT.TFile(fn)
    for j,t in ttree_iterator(f.Get(tn).Get('t')):
        r,l,e,m = t.run, t.lumi, t.event, t.dil_mass
        r = int(r)
        l = int(l)
        e = int(e)
        me[(r,l,e)].append(m)
    return me

def m2s(l):
    return ' * '.join('%.1f' % m for m in l)

def doit2(b1, b2, tn='SimpleNtupler', skip_understood_events=False):
    print 'comparing tree %s in %s to the one in %s' % (tn, b1, b2)
    old = doit('%s/ana_datamc_data.root' % b1, tn)
    new = doit('%s/ana_datamc_data.root' % b2, tn)

    old_json = LumiList('%s/ana_datamc_data.forlumi.json' % b1)
    new_json = LumiList('%s/ana_datamc_data.forlumi.json' % b2)

    oldk = set(old.keys())
    newk = set(new.keys())

    print 'counts:\told\tnew'
    print 'events:\t%i\t%i' % (len(oldk), len(newk))
    print
    print 'diffs in dimu counts?'
    print '%6s%6s%14s%14s%14s%25s%25s' % ('ojs','njs','run','lumi','event','old','new')
    for rle in sorted(oldk & newk):
        if len(old[rle]) != len(new[rle]):
            r,l,e = rle
            print '%6s%6s%14s%14s%14s%25s%25s' % (old_json.contains((r,l)), new_json.contains((r,l)), r,l,e, m2s(old[rle]), m2s(new[rle]))
    print
    print 'symmetric difference:'
    print '%6s%6s%14s%14s%14s%25s%25s' % ('ojs','njs','run','lumi','event','old','new')
    for rle in sorted(oldk.symmetric_difference(newk)):
        r,l,e = rle
        if skip_understood_events and not old_json.contains((r,l)) and new_json.contains((r,l)):
            continue # skip events where it's just because we didn't run on it before
        print '%6s%6s%14s%14s%14s%25s%25s' % (old_json.contains((r,l)), new_json.contains((r,l)), r,l,e, m2s(old[rle]), m2s(new[rle]))

doit2('old/ana_datamc_40ipb', 'ana_datamc_muonsonly')
print '\n************************************************\n'
doit2('old/ana_datamc_34ipb', 'ana_datamc_allgood')
print '\n************************************************\n'
doit2('old/ana_datamc_40ipb', 'ana_datamc_muonsonly', 'SimpleNtuplerSS')
print '\n************************************************\n'
doit2('old/ana_datamc_34ipb', 'ana_datamc_allgood', 'SimpleNtuplerSS')
print '\n************************************************\n'
doit2('old/ana_datamc_34ipb', 'ana_datamc_allgood', 'SimpleNtuplerEmu')
print '\n************************************************\n'


