#!/usr/bin/env python

# (py draw.py >! out.draw) && (py draw.py tight >! out.draw_tight) && tlp plots/nminus1effs*

import sys, os, glob
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)

do_tight = 'tight' in sys.argv
psn = 'plots/nminus1effs'
if do_tight:
    psn += '_tight'
ps = plot_saver(psn, size=(600,600), log=False, pdf=True)

if do_tight:
    nminus1s = [
        #'TiPt',
        'TiDB',
        'TiGlbChi2',
        'TiIso',
        'TiTkHits',
        'TiPxHits',
        'TiMuHits',
        'TiMuMatch',
        ]
else:
    nminus1s = [
        #'NoPt',
        'NoDB',
        'NoGlbChi2',
        'NoIso',
        'NoTkHits',
        'NoPxHits',
        'NoMuHits',
        'NoMuMatch',
        'NoTrgMtch',
        'NoB2B',
        'NoVtxProb',
        'NoCosm',
        ]

mass_ranges = [
    ('60m120',  '60 < m < 120 GeV',  ( 60, 120)),
    ('120m200', '120 < m < 200 GeV', (120, 200)),
    ('200m400', '200 < m < 400 GeV', (200, 400)),
    ('400m',    'm > 400 GeV',       (400, 1e9)),
    ]

to_use = {
    '60m120':  ['ana_nminus1_data.root', 'ana_nminus1_zmumu.root', 'ana_nminus1_ttbar.root'],  #, 'ana_nminus1_inclmu15.root'],
    '120m200': ['ana_nminus1_data.root', 'ana_nminus1_dy120.root', 'ana_nminus1_ttbar.root'],
    '200m400': ['ana_nminus1_data.root', 'ana_nminus1_dy200.root', 'ana_nminus1_ttbar.root'],
    '400m':    ['ana_nminus1_data.root', 'ana_nminus1_dy500.root', 'ana_nminus1_dy1000.root', 'ana_nminus1_ttbar.root'],
    }

pretty = {
    'NoTkHits': '# tk hits > 10',
    'NoPxHits': '# px hits > 0',
    'NoMuStns': '# mu segs > 1',
    'NoDB': '|dB| < 0.2',
    'NoGlbChi2': 'glb #chi^{2}/ndf < 10',
    'NoTkMuon': 'isTrackerMuon',
    'NoMuHits': '# mu hits > 0',
    'NoMuMatch': '# tk. #mu seg > 1',
    'NoCosm': 'anti-cosmic',
    'NoTrgMtch': 'HLT match',
    'NoB2B': 'back-to-back',
    'NoVtxProb': '#chi^{2} #mu#mu vtx < 10',
    'NoIso': 'rel. tk. iso.',
    'ana_nminus1_data.root': 'Data, 366 pb ^{-1}',
    'ana_nminus1_zmumu.root': 'Z#rightarrow#mu#mu, 60 < M < 120 GeV',
    'ana_nminus1_dy120.root': 'DY, M > 120 GeV',
    'ana_nminus1_dy200.root': 'DY, M > 200 GeV',
    'ana_nminus1_dy500.root': 'DY, M > 500 GeV',
    'ana_nminus1_dy1000.root': 'DY, M > 1000 GeV',
    'ana_nminus1_ttbar.root': 't#bar{t}',
    'ana_nminus1_inclmu15.root': 'QCD',
    'ana_nminus1_zssm1000.root': 'Z\' SSM, M=1000 GeV',
    }

styles = {
    'ana_nminus1_data.root':      (ROOT.kBlack,     -1),
    'ana_nminus1_zmumu.root':     (ROOT.kGreen+2, 3001),
    'ana_nminus1_dy120.root':     (ROOT.kGreen+2, 3001),
    'ana_nminus1_dy200.root':     (ROOT.kGreen+2, 3001),
    'ana_nminus1_dy500.root':     (ROOT.kGreen+2, 3001),
    'ana_nminus1_dy1000.root':    (ROOT.kRed,     3002),
    'ana_nminus1_ttbar.root':     (46,            3004),
    'ana_nminus1_inclmu15.root':  (28,            3005),
    }

ymin = {
    '60m120':  0.95,
    '120m200': 0.87,
    '200m400': 0.85,
    '400m':    0.81,
    }
global_ymin = 0.75

def table(f, mass_range):
    hnum = f.Get('NoNo').Get('DileptonMass')
    mlo, mhi = mass_range
    num = get_integral(hnum, mlo, mhi, integral_only=True, include_last_bin=False)
    print 'mass range %5i-%5i:' % mass_range
    print 'numerator:', num
    print '%20s%20s%20s%30s' % ('cut', 'denominator', 'efficiency', '68% CL')
    for nminus1 in nminus1s:
        hden = f.Get(nminus1).Get('DileptonMass')
        den = get_integral(hden, mlo, mhi, integral_only=True, include_last_bin=False)
        e,l,h = clopper_pearson(num, den)
        print '%20s%20f%20f%15f%15f' % (nminus1, den, e, l, h)

ROOT.gStyle.SetTitleX(0.45)
ROOT.gStyle.SetTitleY(0.40)

for name, pretty_name, mass_range in mass_ranges:
    print name, pretty_name

    lg = ROOT.TLegend(0.45, 0.13, 0.94, 0.33)
    lg.SetFillColor(0)
    
    same = 'A'
    effs = []
    
    for fn in to_use[name]:
        print fn
        f = ROOT.TFile(fn)
        table(f, mass_range)

        color, fill = styles[fn]

        l = len(nminus1s)
        nminus1_num = ROOT.TH1F('num', '', l, 0, l)
        nminus1_den = ROOT.TH1F('den', '', l, 0, l)

        hnum = f.Get('NoNo').Get('DileptonMass')
        num = get_integral(hnum, *mass_range, integral_only=True, include_last_bin=False)
        for i,nminus1 in enumerate(nminus1s):
            hden = f.Get(nminus1).Get('DileptonMass')
            den = get_integral(hden, *mass_range, integral_only=True, include_last_bin=False)
            nminus1_num.SetBinContent(i+1, num)
            nminus1_den.SetBinContent(i+1, den)

        eff = binomial_divide(nminus1_num, nminus1_den)
        eff.SetTitle(pretty_name)
        eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin[name], 1.005)
        eff.GetXaxis().SetTitle('cut')
        eff.GetYaxis().SetLabelSize(0.027)
        eff.GetYaxis().SetTitle('n-1 efficiency')
        if 'data' in fn:
            draw = 'P'
            eff.SetLineColor(color)
            lg.AddEntry(eff, pretty.get(fn, fn), 'L')
        else:
            draw = '2'
            eff.SetLineColor(color)
            eff.SetFillColor(color)
            eff.SetFillStyle(fill)
            lg.AddEntry(eff, pretty.get(fn, fn), 'F')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        bnr = eff.GetXaxis().GetNbins()/float(eff.GetN())
        for i in xrange(1,l+1):
            eff.GetXaxis().SetBinLabel(int((i-1)*bnr)+1, pretty.get(nminus1s[i-1], nminus1s[i-1]))
        eff.GetXaxis().LabelsOption('u')

    lg.Draw()
    ps.save(name)
    print
