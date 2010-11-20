#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)
c = ROOT.TCanvas('c', '', 820, 630)
os.system('mkdir -p plots/nminus1effs')

def integ(h,a,b=1e9):
    return h.Integral(h.FindBin(a), h.FindBin(b))

nminus1s = [
    'NoTkHits',
    'NoDB',
    'NoGlbChi2',
    'NoPxHits',
    'NoMuStns',
    'NoTkMuon',
    'NoTrgMtch',
    'NoB2B',
    'NoVtxProb',
]

pretty = {
    'NoTkHits': '# tk hits #geq 10',
    'NoPxHits': '# px hits #geq 1',
    'NoMuStns': '# muon st #geq 2',
    'NoDB': '|d_{xy}| < 0.2',
    'NoGlbChi2': 'glb #chi^{2}/ndf < 10',
    'NoTkMuon': 'isTrackerMuon',
    'NoTrgMtch': 'HLT match',
    'NoB2B': 'B2B cosmics',
    'NoVtxProb': '#chi^{2} #mu#mu vtx < 10',
    'ana_nminus1_data.root': 'Data, 40 pb ^{-1}',
    'ana_nminus1_zmumu.root': 'Z#rightarrow#mu#mu',
    'ana_nminus1_dy120.root': 'Z#rightarrow#mu#mu',
    'ana_nminus1_ttbar.root': 't#bar{t}',
    'ana_nminus1_inclmu15.root': 'QCD',
    'ana_nminus1_zssm750.root': 'Z\' SSM, M=750 GeV',
    }

def table(fn):
    f = ROOT.TFile(fn)
    hnum = f.Get('NoNo').Get('DileptonMass')
    num60120 = integ(hnum, 60, 120)
    num120 = integ(hnum, 120)
    print 'numerator: 60-120:', num60120, ' 120-:', num120
    print '%20s%20s%20s%20s%20s' % ('name', 'den 60-120', 'eff 60-120', 'den 120', 'eff 120')
    for nminus1 in nminus1s:
        hden = f.Get(nminus1).Get('DileptonMass')
        den60120 = integ(hden, 60, 120)
        den120 = integ(hden, 120)
        print '%20s%20f%20f%20f%20f' % (nminus1, den60120, num60120/den60120, den120, num120/den120)

if '120' in sys.argv:
    items = [
        ('ana_nminus1_data.root',    (120, 1e9), ROOT.kBlack, -1),
        ('ana_nminus1_dy120.root',   (120, 1e9), ROOT.kBlue, 3006),
        ('ana_nminus1_ttbar.root',   (120, 1e9), ROOT.kRed, 3007),
        ('ana_nminus1_zssm750.root', (120, 1e9), ROOT.kGreen+2, 3004)
        ]
    lg = ROOT.TLegend(0.14, 0.14, 0.63, 0.34)
else:         
    items = [
        ('ana_nminus1_data.root',     (60,120), ROOT.kBlack, -1),
        ('ana_nminus1_zmumu.root',    (60,120), ROOT.kBlue, 3006),
        ('ana_nminus1_ttbar.root',    (60,120), ROOT.kRed, 3007),
        ]
    lg = ROOT.TLegend(0.14, 0.14, 0.50, 0.32)

same = 'A'
effs = []
lg.SetFillColor(0)
for fn, mass_range, color, fill in items:
    f = ROOT.TFile(fn)

    l = len(nminus1s)
    nminus1_num = ROOT.TH1F('num', '', l, 0, l)
    nminus1_den = ROOT.TH1F('den', '', l, 0, l)

    hnum = f.Get('NoNo').Get('DileptonMass')
    num = integ(hnum, *mass_range)
    for i,nminus1 in enumerate(nminus1s):
        hden = f.Get(nminus1).Get('DileptonMass')
        den = integ(hden, *mass_range)
        nminus1_num.SetBinContent(i+1, num)
        nminus1_den.SetBinContent(i+1, den)

    eff = binomial_divide(nminus1_num, nminus1_den)
    if '120' in sys.argv:
        eff.SetTitle('M_{#mu#mu} > 120 GeV')
    else:
        eff.SetTitle('60 < M_{#mu#mu} < 120 GeV')
    eff.GetXaxis().SetRangeUser(0,12)
    eff.GetYaxis().SetRangeUser(0.95,1.005)
    eff.GetXaxis().SetTitle('cut')
    eff.GetYaxis().SetLabelSize(0.03)
    eff.GetYaxis().SetTitle('n-1 efficiency')
    eff.GetYaxis().SetTitle('n-1 efficiency')
    if 'data' in fn:
        draw = 'P'
        eff.SetLineColor(color)
        lg.AddEntry(eff, pretty[fn], 'L')
    else:
        draw = '2'
        eff.SetLineColor(color)
        eff.SetFillColor(color)
        eff.SetFillStyle(fill)
        lg.AddEntry(eff, pretty[fn], 'F')
    draw += same
    eff.Draw(draw)
    effs.append(eff)
    same = ' same'
    bnr = eff.GetXaxis().GetNbins()/eff.GetN()
    for i in xrange(1,l+1):
        eff.GetXaxis().SetBinLabel((i-1)*bnr+1, pretty[nminus1s[i-1]])
    eff.GetXaxis().LabelsOption('u')

    print fn
    table(fn)

lg.Draw()
fn = 'plots/nminus1effs/nminus1%s.png' % ('120' if '120' in sys.argv else '')
c.SaveAs(fn)
c.SaveAs(fn.replace('png', 'root'))

