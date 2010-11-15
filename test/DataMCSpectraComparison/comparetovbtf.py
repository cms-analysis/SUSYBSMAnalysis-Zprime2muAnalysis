#!/usr/bin/env python

from collections import defaultdict
from math import sinh, cos, sin
from FWCore.PythonUtilities.LumiList import LumiList
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, ttree_iterator, set_zp2mu_style
set_zp2mu_style()

me = defaultdict(list)
vbtf = defaultdict(list)

me_json = LumiList('ana_datamc/ana_datamc_data.forlumi.json')
vbtf_json = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON.txt')

peak_only = False
tracker_mass = True

def ptetaphi2p4(pt,eta,phi):
    px = pt * cos(phi)
    py = pt * sin(phi)
    pz = pt * sinh(eta)
    e = (pt**2 + pz**2 + 0.10566**2)**0.5
    return px,py,pz,e

fns = ['ana_datamc/ana_datamc_data.root']
for fn in fns:
    f = ROOT.TFile(fn)
    for j,t in ttree_iterator(f.SimpleNtuplerVBTF.Get('t')):
        if not peak_only or (t.dil_mass > 60 and t.dil_mass < 120):
            r,l,e,m = t.run, t.lumi, t.event, t.dil_mass
            if tracker_mass:
                px0,py0,pz0,e0 = ptetaphi2p4(t.lep_tk_pt[0], t.lep_tk_eta[0], t.lep_tk_phi[0])
                px1,py1,pz1,e1 = ptetaphi2p4(t.lep_tk_pt[1], t.lep_tk_eta[1], t.lep_tk_phi[1])
                m = ((e0+e1)**2 - (px0+px1)**2 - (py0+py1)**2 - (pz0+pz1)**2)**0.5
            if m < 1:
                print 'mass < 1 GeV:', m, t.lep_tk_pt[0], t.lep_tk_eta[0], t.lep_tk_phi[0], t.lep_tk_pt[1], t.lep_tk_eta[1], t.lep_tk_phi[1]
            if m < 20 and m > 1:
                print '    mass < 20 GeV:', m, t.lep_tk_pt[0], t.lep_tk_eta[0], t.lep_tk_phi[0], t.lep_tk_pt[1], t.lep_tk_eta[1], t.lep_tk_phi[1]
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

mass_diffs = {}
h_mass_diff = ROOT.TH1F('h_mass_diff', '', 10000, -30, 30)

def m2s(l):
    return ' * '.join('%.1f' % m for m in l)

print 'counts:\tme\tvbtf'
print 'events:\t%i\t%i' % (len(mek), len(vbtfk))
print
print 'diffs in dimu counts?'
print '%14s%14s%14s%25s%25s' % ('run','lumi','event', 'me', 'vbtf')
for rle in sorted(mek & vbtfk):
    if len(me[rle]) != len(vbtf[rle]):
        print '%14s%14s%14s%25s%25s' % (rle + (m2s(me[rle]), m2s(vbtf[rle])))
    else:
        assert(len(me[rle]) == 1)
        mass_diffs[rle] = dm = me[rle][0] - vbtf[rle][0]
        h_mass_diff.Fill(dm)
print
h_mass_diff.Draw()
ROOT.c1.SetLogy(1)
ROOT.c1.SaveAs('~/asdf/massdiffs.png')
ROOT.c1.SaveAs('~/asdf/massdiffs.root')
print 'max mass diff:', max(abs(m) for m in mass_diffs.values())
print
print '|delta m| > 0.01 GeV:\n'
print '%14s%14s%14s%25s%25s' % ('run','lumi','event', 'me', 'vbtf')
for rle in sorted(mass_diffs.keys()):
    if mass_diffs[rle] > 0.01:
        print '%14s%14s%14s%25s%25s' % (rle + (m2s(me[rle]), m2s(vbtf[rle])))
print
print 'symmetric difference where rle is in both jsons:'
print '%14s%14s%14s%25s%25s' % ('run','lumi','event', 'me', 'vbtf')
only_in_me_json = defaultdict(list)
only_in_vbtf_json = defaultdict(list)
for rle in sorted(mek.symmetric_difference(vbtfk)):
    r,l,e = rle
    if not me_json.contains((r,l)):
        only_in_vbtf_json[rle].append(vbtf[rle])
        continue
    if not vbtf_json.contains((r,l)):
        only_in_me_json[rle].append(me[rle])
        continue
    print '%14s%14s%14s%25s%25s' % (r,l,e, m2s(me[rle]), m2s(vbtf[rle]))
print
print 'events only in my json:', len(only_in_me_json)
print 'events only in vbtf json:', len(only_in_vbtf_json)
