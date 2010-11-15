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
    
def f(which_plot, xtitle, ytitle, m=1000):
    f = ROOT.TFile('ana_effres_zp%i.root' % m)
    cache = []
    tracks = ['inner', 'tpfms', 'picky', 'pmc', 'tmr', 'sigmaswitch']
    histos = {}
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
        #h.SetStats(0)
        fit_gaussian(h, 1.5, draw=True)
        tl = ROOT.TLatex(0.13, 0.90, nice_name[t])
        tl.SetTextSize(0.06)
        tl.SetNDC(1)
        cache.append(tl)
        tl.Draw()
        p.Update()
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

f('LeptonInvPtRes', 'q/p_{T} relative resolution', 'Entries/0.01')
f('LeptonPtRes', 'q*p_{T} relative resolution', 'Entries/0.01')
f('DileptonMassRes', 'm_{#mu#mu} relative resolution', 'Entries/0.01')
