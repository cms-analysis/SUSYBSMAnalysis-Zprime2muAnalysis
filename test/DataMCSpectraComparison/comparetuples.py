#!/usr/bin/env python

from collections import defaultdict
from FWCore.PythonUtilities.LumiList import LumiList
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, detree

def doit(fn, mass_name, dir_name, tree_name, cut):
    f = ROOT.TFile(fn)
    t = f.Get(dir_name).Get(tree_name)
    rlem = detree(t, 'run:lumi:event:' + mass_name, cut, lambda x: (int(x[0]), int(x[1]), int(x[2]), float(x[3])))
    me = defaultdict(list)
    for r,l,e,m in rlem:
        me[(r,l,e)].append(m)
    return me

def m2s(l):
    if not hasattr(l, '__iter__'):
        l = [l]
    return ' * '.join('%.1f' % m for m in l)


def doit2(skip_understood_events=False, max_mass_diff=0.1, use_jsons=True):
    old = doit('data/Run2011MuonsOnly/ana_datamc_data.root', 'dil_mass', 'SimpleNtupler', 't', 'OurSelNewNoSign && SameSign')
    new = doit('m10t/ntuple_SingleMuRun2011.root',           'mass',     'SameSign',       't', '')
    #old = doit('mc/ana_datamc_ttbar.root',  'dil_mass', 'SimpleNtupler', 't', 'OurSelNewNoSign && SameSign')
    #new = doit('m10t/ntuple_MC_ttbar.root', 'mass',     'SameSign',       't', '')
    #use_jsons = False

    if use_jsons:
        old_json = LumiList('data/Run2011MuonsOnly/ana_datamc_data.forlumi.json')
        new_json = LumiList('m10t/ntuple_SingleMuRun2011.report.json')
    else:
        class dummy:
            def contains(*args):
                return 'N/A'
        old_json = new_json = dummy()

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
        mold = max(old[rle] + [0])
        mnew = max(new[rle] + [0])
        #if not (mold > 300 or mnew > 300):
        #    continue
        if skip_understood_events and not old_json.contains((r,l)) and new_json.contains((r,l)):
            continue # skip events where it's just because we didn't run on it before
        print '%6s%6s%14s%14s%14s%25s%25s' % (old_json.contains((r,l)), new_json.contains((r,l)), r,l,e, m2s(old[rle]), m2s(new[rle]))

doit2()
