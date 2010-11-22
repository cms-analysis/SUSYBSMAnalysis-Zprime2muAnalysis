#!/usr/bin/env python

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
    'pmc': 'Tune P',
    'tmr': 'TMR',
    'sigmaswitch': 'Sigma-switch'
    }
    
def f(which_plot, xtitle, ytitle, m=1000, no_stats=True, errors_in_table=False):
    print which_plot
    f = ROOT.TFile('ana_effres_zp%i.root' % m)
    cache = []
    tracks = ['inner', 'tpfms', 'picky', 'pmc', 'tmr', 'sigmaswitch']
    histos = {}
    print '%12s%10s%10s%10s%25s%25s%25s%25s%25s' % ('reco', 'entries', 'under', 'over', 'mean', 'rms', 'constant', 'mu', 'sigma')
    for i,t in enumerate(tracks):
        d = f.Get('Resolution%s' % t)
        h = d.Get(which_plot)
        histos[t] = h
        p = c.cd(i+1)
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.GetYaxis().SetLabelSize(0.02)
        #h.SetTitle(t)
        h.SetTitle('')
        if no_stats:
            h.SetStats(0)
        stats = get_hist_stats(h, factor=1.5, draw=True)
        stats['reco'] = t
        for x in ['mean', 'rms', 'constant', 'mu', 'sigma']:
            if errors_in_table:
                stats[x] = '%6.4f \pm %6.4f' % stats[x]
            else:
                stats[x] = '%10.4f' % stats[x][0]
        print '%(reco)12s%(entries)10i%(under)10i%(over)10i%(mean)25s%(rms)25s%(constant)25s%(mu)25s%(sigma)25s' % stats
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
