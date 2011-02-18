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
    if not hasattr(l, '__iter__'):
        l = [l]
    return ' * '.join('%.1f' % m for m in l)

def doit2(b1, b2, tn='SimpleNtupler', skip_understood_events=False, max_mass_diff=0.1):
    print 'comparing tree %s in %s to the one in %s' % (tn, b1, b2)
    old = doit('%s/ana_datamc_data.root' % b1, tn)
    new = doit('%s/ana_datamc_data.root' % b2, tn)

    old_json = LumiList('%s/ana_datamc_data.forlumi.json' % b1)
    new_json = LumiList('%s/ana_datamc_data.forlumi.json' % b2)

    oldk = set(old.keys())
    newk = set(new.keys())
    max_mass_diff_seen = 0
    
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
    print 'mass comparison:'
    print '%6s%6s%14s%14s%14s%25s%25s%10s' % ('ojs','njs','run','lumi','event','old','new','maxd')
    to_print = []
    for rle in sorted(oldk & newk):
        if len(old[rle]) == len(new[rle]):
            if len(old[rle]) > 1:
                print '***\nwarning: more than one dilepton in mass comparison, just dumbly sorting\n***'
            oldm = sorted(old[rle])
            newm = sorted(new[rle])
            masses_diff = [abs(o-n) for o,n in zip(oldm,newm)]
            mmd = max(masses_diff)
            if mmd > max_mass_diff_seen:
                max_mass_diff_seen = mmd
            if mmd > max_mass_diff:
                r,l,e = rle
                to_print.append((mmd, '%6s%6s%14s%14s%14s%25s%25s%10s' % (old_json.contains((r,l)), new_json.contains((r,l)), r,l,e, m2s(old[rle]), m2s(new[rle]), m2s(mmd))))
    to_print.sort(key=lambda x: x[0], reverse=True)
    for p in to_print:
        print p[1]
    print 'max mass diff seen: %.5f' % max_mass_diff_seen
    print
    print 'symmetric difference:'
    print '%6s%6s%14s%14s%14s%25s%25s' % ('ojs','njs','run','lumi','event','old','new')
    for rle in sorted(oldk.symmetric_difference(newk)):
        r,l,e = rle
        if skip_understood_events and not old_json.contains((r,l)) and new_json.contains((r,l)):
            continue # skip events where it's just because we didn't run on it before
        print '%6s%6s%14s%14s%14s%25s%25s' % (old_json.contains((r,l)), new_json.contains((r,l)), r,l,e, m2s(old[rle]), m2s(new[rle]))

if 'simple' in sys.argv:
    n = sys.argv.index('simple')
    doit2(sys.argv[n+1], sys.argv[n+2])
    raise 'done'

doit2('ana_datamc_current/muonsonly', 'ana_datamc_nov4/muonsonly')
print '\n************************************************\n'
doit2('ana_datamc_current/muonsonly', 'ana_datamc_nov4/muonsonly', 'SimpleNtuplerSS')
print '\n************************************************\n'
doit2('ana_datamc_current/muonsonly', 'ana_datamc_nov4/muonsonly', 'SimpleNtuplerVBTF')
print '\n************************************************\n'
raise 'done'

doit2('old/ana_datamc_34ipb', 'ana_datamc_allgood')
print '\n************************************************\n'
doit2('old/ana_datamc_40ipb', 'ana_datamc_muonsonly', 'SimpleNtuplerSS')
print '\n************************************************\n'
doit2('old/ana_datamc_34ipb', 'ana_datamc_allgood', 'SimpleNtuplerSS')
print '\n************************************************\n'
doit2('old/ana_datamc_34ipb', 'ana_datamc_allgood', 'SimpleNtuplerEmu')
print '\n************************************************\n'


