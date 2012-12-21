#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

from pprint import pprint
import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)
ROOT.TH1.AddDirectory(0)

outfile = ROOT.TFile("whargl.root","recreate")
iarp=0
do_tight = 'tight' in sys.argv
psn = 'plots/nminus1effs'
if do_tight:
    psn += '_tight'
ps = plot_saver(psn, size=(600,600), log=False, pdf=True)
ps.c.SetBottomMargin(0.2)

if do_tight:
    nminus1s = [
        #'TiPt',
        'TiDB',
        'TiGlbChi2',
        'TiIso',
        'TiTkLayers',
        'TiPxHits',
        'TiMuHits',
        'TiMuMatch',
        ]
else:
    nminus1s = [
        'NoPt',
        'NoDB',
        'NoIso',
        'NoTkLayers',
        'NoPxHits',
        'NoMuHits',
        'NoMuMatch',
        'NoVtxProb',
        'NoB2B',
        'NoDptPt',
        'NoCosm',
        'NoTrgMtch',
        ]

pretty = {
    'NoPt': 'p_{T} > 45 GeV',
    'NoTkLayers': '# tk lay > 5',
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
    'NoDptPt': 'dpT/pT',
    'NoIso': 'rel. tk. iso.',
    'data': 'Data, %.1f fb^{-1}',
    'mcsum': 'Simulation',
    'zmumu': 'Z#rightarrow#mu#mu, 60 < M < 120 GeV',
    'dy120_c1': 'DY#rightarrow#mu#mu, M > 120 GeV',
    'dy200_c1': 'DY#rightarrow#mu#mu, M > 200 GeV',
    'dy500_c1': 'DY#rightarrow#mu#mu, M > 500 GeV',
    'dy1000_c1': 'DY#rightarrow#mu#mu, M > 1000 GeV',
    'ttbar': 't#bar{t}',
    'inclmu15': 'QCD',
    'zssm1000': 'Z\' SSM, M=1000 GeV',
    '60m120': '60 < m < 120 GeV',  
    '120m200': '120 < m < 200 GeV', 
    '200m400': '200 < m < 400 GeV',
    '400m600': '400 < m < 600 GeV',
    '600m': 'm > 600 GeV',
    }

class nm1entry:
    def __init__(self, sample, is_data=False):
        if type(sample) == str:
            self.name = sample
            self.fn = self.make_fn(sample) if is_data else None
        else:
            self.name = sample.name
            self.fn = self.make_fn(self.name)
        self.prepare_histos()
            
    def make_fn(self, name):
        return 'nminus1_histos/ana_nminus1_%s.root' % name
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                self.histos[nminus1] = f.Get(nminus1).Get('DileptonMass').Clone()

    def prepare_histos_sum(self, samples, lumi):
        self.histos = {}
        for nminus1 in nminus1s + ['NoNo']:
            hs = []
            for sample in samples:
                f = ROOT.TFile(self.make_fn(sample.name))
                h = f.Get(nminus1).Get('DileptonMass').Clone()
                print nminus1, sample.name, sample.partial_weight*lumi
                h.Scale(sample.partial_weight * lumi)
                hs.append(h)
            hsum = hs[0].Clone()
            for h in hs[1:]:
                hsum.Add(h)
            self.histos[nminus1] = hsum

data, lumi = nm1entry('data', True), 15796. # lumi in pb
mcsum = nm1entry('mcsum')

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import zmumu, dy120_c1, dy200_c1, dy500_c1, dy1000_c1, ttbar, inclmu15
raw_samples = [zmumu, dy120_c1, dy200_c1, dy500_c1, dy1000_c1, ttbar, inclmu15]
mcsum.prepare_histos_sum(raw_samples, lumi)
mc_samples = [nm1entry(sample) for sample in raw_samples]
for mc_sample in mc_samples:
    exec '%s = mc_sample' % mc_sample.name

mass_ranges = [
    ('60m120',  ( 60, 120)),
    ('120m200', (120, 200)),
    ('200m400', (200, 400)),
    ('400m600', (400, 600)),
    ('600m',    (600, 1e9)),
    ]

to_use = {
    '60m120':  [data, mcsum], #, zmumu, ttbar],
    '120m200': [data, mcsum], #, dy120_c1, ttbar],
    '200m400': [data, mcsum], #, dy200_c1, ttbar],
    '400m600': [data, mcsum], #, dy200_c1, ttbar],
    '600m':    [data, mcsum], #, dy500_c1, ttbar],
    }

styles = {
    'data':      (ROOT.kBlack,     -1),
    'mcsum':     (ROOT.kMagenta, 1001),
    'zmumu':     (ROOT.kRed,     1001),
    'dy120_c1':  (ROOT.kRed,     1001),
    'dy200_c1':  (ROOT.kRed,     1001),
    'dy500_c1':  (ROOT.kRed,     1001),
    'dy1000_c1': (ROOT.kBlue,    1001),
    'ttbar': (ROOT.kGreen+2, 3005),
    'inclmu15':  (28,            3002),
    }

ymin = {
    '60m120':  0.95,
    '120m200': 0.87,
    '200m400': 0.85,
    '600m':    0.81,
    }
global_ymin = 0.58

def table(entry, mass_range):
    print entry.name
    hnum = entry.histos['NoNo']
    mlo, mhi = mass_range
    num = get_integral(hnum, mlo, mhi, integral_only=True, include_last_bin=False)
    print 'mass range %5i-%5i:' % mass_range
    print 'numerator:', num
    print '%20s%20s%21s%21s' % ('cut', 'denominator', 'efficiency', '68% CL')
    for nminus1 in nminus1s:
        hden = entry.histos[nminus1]
        den = get_integral(hden, mlo, mhi, integral_only=True, include_last_bin=False)
        e,l,h = clopper_pearson(num, den)
        print '%20s%20f%20f%15f%15f' % (nminus1, den, e, l, h)

ROOT.gStyle.SetTitleX(0.25)
ROOT.gStyle.SetTitleY(0.50)

for name, mass_range in mass_ranges:
    pretty_name = pretty[name]
    print name, pretty_name

    lg = ROOT.TLegend(0.25, 0.23, 0.81, 0.43)
    lg.SetTextSize(0.05)
    lg.SetFillColor(0)
    lg.SetBorderSize(1)
    
    same = 'A'
    effs = []
    
    for entry in to_use[name]:
        table(entry, mass_range)
        color, fill = styles[entry.name]

        l = len(nminus1s)
        nminus1_num = ROOT.TH1F('num', '', l, 0, l)
        nminus1_den = ROOT.TH1F('den', '', l, 0, l)

        hnum = entry.histos['NoNo']
        num = get_integral(hnum, *mass_range, integral_only=True, include_last_bin=False)
        for i,nminus1 in enumerate(nminus1s):
            hden = entry.histos[nminus1]
            den = get_integral(hden, *mass_range, integral_only=True, include_last_bin=False)
            nminus1_num.SetBinContent(i+1, num)
            nminus1_den.SetBinContent(i+1, den)

        eff = binomial_divide(nminus1_num, nminus1_den)
        eff.SetTitle(pretty_name)
        eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin[name], 1.01)
        eff.GetXaxis().SetTitle('cut')
        eff.GetYaxis().SetLabelSize(0.027)
        eff.GetYaxis().SetTitle('n-1 efficiency')
        if 'data' in entry.name:
            draw = 'P'
            eff.SetLineColor(color)
            eff.SetMarkerStyle(20)
            eff.SetMarkerSize(1.05)
            lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
        else:
            draw = '2'
            eff.SetLineColor(color)
            eff.SetFillColor(color)
            eff.SetFillStyle(fill)
            lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'L')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        bnr = eff.GetXaxis().GetNbins()/float(eff.GetN()+1)
        for i, n in enumerate(nminus1s):
            eff.GetXaxis().SetBinLabel(int((i+0.5)*bnr), pretty.get(n,n))
        eff.GetXaxis().LabelsOption('v')
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1

    lg.Draw()
    ps.save(name)
    print
