#!/usr/bin/env python

# (py diffrecores.py >! plots/out.diffrecores) && tlp plots/*diffrecores

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.0)

os.system('mkdir -p plots/diffrecores')
c = ROOT.TCanvas('c', '', 820*2, 630*2)
c.Divide(3,2)

nice_name = {
    'inner': 'Tracker-only',
    'tpfms': 'TPFMS',
    'picky': 'Picky',
    'tunep': 'Tune P',
    'tmr': 'TMR',
    'global': 'Global',
    'sigmaswitch': 'Sigma-switch',
    'LeptonInvPtRes': r'Muon $q / \pt$',
    'LeptonPtRes': r'Muon $q \cdot \pt$',
    'DileptonMassRes': 'Dimuon invariant mass',
    }
    
def f(which_plot, xtitle, ytitle, m=1000, no_stats=True, errors_in_table=True):
    print r'\hline\hline'
    print r'\multicolumn{6}{|l|}{%s} \\' % nice_name[which_plot]
    f = ROOT.TFile('ana_effres_zp%i.root' % m)
    cache = []
    tracks = ['inner', 'tpfms', 'picky', 'tunep', 'tmr', 'sigmaswitch', 'global']
    histos = {}
    print r'\hline'
    print r'%20s & %15s & %25s & %25s & %15s & %15s \\' % ('Fit/selector', 'Entries', r'Fitted $\sigma$, \%', r'RMS, \%', r'$\# < -0.5$', r'$\# > 0.5$')
    print r'\hline'
    for i,t in enumerate(tracks):
        draw = t != 'global'
        d = f.Get('Resolution%s' % t)
        h = d.Get(which_plot)
        histos[t] = h
        if draw:
            p = c.cd(i+1)
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.GetYaxis().SetLabelSize(0.02)
        #h.SetTitle(t)
        h.SetTitle('')
        if no_stats:
            h.SetStats(0)
        stats = get_hist_stats(h, factor=1.5, draw=draw)
        stats['reco'] = nice_name[t]
        for x in ['mean', 'rms', 'constant', 'mu', 'sigma']:
            if errors_in_table:
                y,z = stats[x]
                y *= 100
                z *= 100
                stats[x] = '%6.1f $\pm$ %6.2f' % (y,z)
            else:
                stats[x] = '%6.1f' % stats[x][0]*100
        print '%(reco)20s & %(entries)15i & %(sigma)25s & %(rms)25s & %(under)15i & %(over)15i \\\\' % stats
        if not draw:
            continue
        tl = ROOT.TLatex(0.13, 0.90, nice_name[t])
        tl.SetTextSize(0.06)
        tl.SetNDC(1)
        cache.append(tl)
        tl.Draw()
        p.Update()
        if not no_stats:
            s = h.FindObject('stats')
            s.SetX1NDC(0.39)
            s.SetY1NDC(0.37)
            s.SetX2NDC(0.99)
            s.SetY2NDC(0.77)
            s.SetOptStat(111110)
            s.SetOptFit(11)
            s.Draw()
    c.SaveAs('plots/diffrecores/%s%i.png' % (which_plot, m))
    c.SaveAs('plots/diffrecores/%s%i.root' % (which_plot, m))
    print

f('LeptonInvPtRes', 'q/p_{T} relative resolution', 'Entries/0.01')
f('LeptonPtRes', 'q*p_{T} relative resolution', 'Entries/0.01')
f('DileptonMassRes', 'm_{#mu#mu} relative resolution', 'Entries/0.01')
