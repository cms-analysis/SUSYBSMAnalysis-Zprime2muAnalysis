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
#psn = 'plots/nminus1effs'
psn = 'plots/paper2016'
# 'plots' = '/afs/cern.ch/work/c/cschnaib/NMinus1Effs/plots/TAG/'
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
                #'NoCosm',
        'NoTrgMtch',
        ]

pretty = {
    'NoPt': 'p_{T} > 53 GeV',
    'NoTkLayers': '# tk lay > 5',
    'NoPxHits': '# px hits > 0',
    'NoMuStns': '# mu segs > 1',
    'NoDB': '|dxy| < 0.2',
    'NoGlbChi2': 'glb #chi^{2}/ndf < 10',
    'NoTkMuon': 'isTrackerMuon',
    'NoMuHits': '# mu hits > 0',
    'NoMuMatch': '# matched stations > 1',
    'NoCosm': 'anti-cosmic',
    'NoTrgMtch': 'HLT match',
    'NoB2B': 'back-to-back',
    'NoVtxProb': '#chi^{2} #mu#mu vtx < 20',
    'NoDptPt': 'dpT/pT',
    'NoIso': 'rel. tk. iso.',
    #'data': 'Data, %.1f fb^{-1}',
    'data': 'Data, %.1f fb^{-1}, MuonOnly',
    'dataB': 'Data RunB, %.1f fb^{-1}, MuonOnly',
    'dataCD': 'Data RunC+D, %.1f fb^{-1}, MuonOnly',
    'dataBCD': 'Data RunB+C+D, %.1f fb^{-1}, MuonOnly',
    'mcsum_lumi': 'Simulation',
    'mcsum_ref': 'Simulation',
    'mc50m120_lumi': 'Simulation 60 < M < 120 GeV',
    'mc50m120_ref': 'Simulation 60 < M < 120 GeV',
    'mc120m_lumi': 'Simulation M > 120 GeV',
    'mc120m_ref': 'Simulation M > 120 GeV',
    'mc800m2300_lumi': 'Simulation 800 < M < 2300 GeV',
    'mc800m2300_ref': 'Simulation 800 < M < 2300 GeV',
    'mc400m2300_lumi': 'Simulation 400 < M < 2300 GeV',
    'mc400m2300_ref': 'Simulation 400 < M < 2300 GeV',
    'zmumu': 'Z#rightarrow#mu#mu, 60 < M < 120 GeV',
    'dy120_c1': 'DY#rightarrow#mu#mu, M > 120 GeV',
    'dy200_c1': 'DY#rightarrow#mu#mu, M > 200 GeV',
    'dy500_c1': 'DY#rightarrow#mu#mu, M > 500 GeV',
    'dy1000_c1': 'DY#rightarrow#mu#mu, M > 1000 GeV',
    'dy50': 'DY#rightarrow#mu#mu madgraph',
#    'dy50': 'DY#rightarrow#mu#mu, M > 50 GeV',
    'dy50to120': 'DY#rightarrow#mu#mu powheg',
    'dy50_startup': 'DY#rightarrow#mu#mu startup',
    'ttbar': 't#bar{t}',
    'ttbar_pow': 't#bar{t} powheg',
    'ttbar_startup': 't#bar{t} startup',
    'ww_incl': 'WW',
    'zz_incl': 'ZZ',
    'wz' : 'WZ',
    'inclmu15': 'QCD',
    'zssm1000': 'Z\' SSM, M=1000 GeV',
    'zpsi5000': 'Z\'_{#psi}, M=5000 GeV',
    'zpsi5000_m1TeV': 'Z\'_{#psi}, M=5000 GeV',
    'zpsi5000_1m3TeV': 'Z\'_{#psi}, M=5000 GeV',
    'zpsi5000_3mTeV': 'Z\'_{#psi}, M=5000 GeV',
    '60m120_BCD': '60 < m < 120 GeV',
    '60m120_CD': '60 < m < 120 GeV',
    '60m120': '60 < m < 120 GeV',
    '60m': 'm > 60 GeV',
    '70m110': '70 < m < 110 GeV',
    '120m200': '120 < m < 200 GeV', 
    '200m400': '200 < m < 400 GeV',
    '400m600': '400 < m < 600 GeV',
    '200m': 'm > 200 GeV',
    '50m': 'm > 50 GeV',
    '70m': 'm > 70 GeV',
    '120m_BCD': 'm > 120 GeV',
    '120m_CD': 'm > 120 GeV',
    '120m': 'm > 120 GeV',
    'DY120to200Powheg': 'DY#rightarrow#mu#mu 120 < m < 200 GeV',
    'DY200to400Powheg': 'DY#rightarrow#mu#mu 200 < m < 400 GeV',
    'DY400to800Powheg': 'DY#rightarrow#mu#mu 400 < m < 800 GeV',
    'DY800to1400Powheg': 'DY#rightarrow#mu#mu 800 < m < 1400 GeV',
    'dy1400to2300': 'DY#rightarrow#mu#mu 1400 < m < 2300 GeV',
    '400m800' : '400 < m < 800 GeV',
    '800m1400': '800 < m < 1400 GeV',
    '1400m2300':'1400 < m < 2300 GeV',
    '800m2300':'800 < m < 2300 GeV',
    '400m2300':'400 < m < 2300 GeV',
    'all_lumi':'Simulation M > 120 GeV',
    'all_ref':'Simulation M > 120 GeV',
    '120m1400':'120 < M < 1400 GeV',
    }

class nm1entry:
    def __init__(self, sample, is_data, lumi):
        if type(sample) == str:
            self.name = sample
            self.fn = 'data/ana_nminus1_%s.root' %sample
            #self.fn = self.make_fn(sample) if is_data else None
            self.lumi = lumi if is_data else None
            self.is_data = is_data
        else:
            self.name = sample.name
            self.fn = self.make_fn(self.name)
            self.partial_weight = sample.partial_weight
            self.is_data = is_data
        self.prepare_histos()
            

    def make_fn(self, name):
        '''
        if self.is_data:
            return 'data/ana_nminus1_%s.root' % name
        else:
        '''
        return 'mc/ana_nminus1_%s.root' % name
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                if nminus1=='NoVtxProb':
                    self.histos[nminus1] = f.Get(nminus1).Get('DileptonMass').Clone()
                else:
                    self.histos[nminus1] = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()#DileptonMass

    def prepare_histos_sum(self, samples, lumi):
        self.histos = {}
        for nminus1 in nminus1s + ['NoNo']:
            hs = []
            #print '%20s%20s%21s%20s%20s' % ('cut', 'sampe name', 'partial weight', 'scale(ref)','scale(lumi)')
            for sample in samples:
                f = ROOT.TFile(self.make_fn(sample.name))
                if nminus1 == 'NoVtxProb':
                    h = f.Get(nminus1).Get('DileptonMass').Clone()
                else:
                    h = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()
                #print '%20s%20s%20.15f%20f%20f' % (nminus1, sample.name, sample.partial_weight, refN/refXS, lumiBCD)
                # partial_weight = cross_section * k_factor / Nevents
                if lumi>0:
                    # scale to luminosity for comparision of single dataset to MC
                    h.Scale(sample.partial_weight * lumi) 
                    #print nminus1, sample.name, sample.partial_weight*lumi
                    #print '%20s%20s%20.10f%20f' % (nminus1, sample.name, sample.partial_weight, lumi)
                if lumi<0:
                    # scale to reference cross section/Nevents for comparision of multiple datasets to MC
                    h.Scale(sample.partial_weight * refN / refXS)  
                    #print nminus1, sample.name, sample.partial_weight*refN/refXS
                    #print '%20s%20s%20.10f%20f' % (nminus1, sample.name, sample.partial_weight, refN/refXS)
                hs.append(h)
            hsum = hs[0].Clone()
            for h in hs[1:]:
                hsum.Add(h)
            self.histos[nminus1] = hsum

#data, lumi = nm1entry('data', True), 242.8 # lumi in pb
nolumi = -1
lumiB = 50.7 
lumiCD = 2619.44
lumiD = 2572.19
#lumiBCD = 2660.14
lumiBCD = 2800.
#dataB = nm1entry('dataB', True, lumiB)#lumiB )
#dataCD = nm1entry('dataCD', True, lumiCD)#lumiCD )
dataBCD = nm1entry('dataBCD', True, lumiBCD)#lumiCD )
#dataD = nm1entry('dataD', True, lumiD )
#mcsum_lumi = nm1entry('mcsum_lumi',False,lumiBCD)
#mcsum_ref = nm1entry('mcsum_ref',False,nolumi)
#mc50m_lumi = nm1entry('mc50m_lumi',False,lumiBCD)
#mc50m_ref = nm1entry('mc50m_ref',False,nolumi)
#mc50m120_lumi = nm1entry('mc50m120_lumi',False,lumiBCD)
#mc50m120_ref = nm1entry('mc50m120_ref',False,nolumi)
#mc120m_lumi = nm1entry('mc120m_lumi',False,lumiBCD)
#mc120m_ref = nm1entry('mc120m_ref',False,nolumi)
#mc800m2300_lumi = nm1entry('mc800m2300_lumi',False,lumiCD)
#mc800m2300_ref = nm1entry('mc800m2300_ref',False,nolumi)
#mc400m2300_lumi = nm1entry('mc400m2300_lumi',False,lumiCD)
#mc400m2300_ref = nm1entry('mc400m2300_ref',False,nolumi)

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
# old startup alignment
#raw_samples = [dy50to120,DY120to200Powheg,DY200to400Powheg,DY400to800Powheg,DY800to1400Powheg,dy1400to2300,dy2300to3500,DY3500to4500Powheg,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,wjets,inclmu15,zpsi5000,qcd50to80,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]
#raw_samples = [dy50to120,DY120to200Powheg,DY200to400Powheg,DY400to800Powheg,DY800to1400Powheg,dy1400to2300,dy2300to3500,DY3500to4500Powheg,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,zpsi5000]
#raw_samples50m120 = [dy50to120,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,wjets,qcd50to80,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]#,inclmu15]
#raw_samples120m = [DY120to200Powheg,DY200to400Powheg,DY400to800Powheg,DY800to1400Powheg,dy2300to3500,DY3500to4500Powheg,dy4500to6000,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,wjets,qcd50to80,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]#inclmu15
#qcd120to170,qcd600to800
# new startup alignment

#raw_samples = [dy50to120_s,dy120to200_s,dy200to400_s,dy400to800_s,dy800to1400_s,dy1400to2300_s,dy2300to3500_s,dy3500to4500_s,dy4500to6000_s,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,wjets,inclmu15,zpsi5000,qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]#,dy6000_s
#raw_samples50m120 = [dy50to120_s,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]#,inclmu15],dy6000_s
#raw_samples120m = [dy120to200_s,dy200to400_s,dy400to800_s,dy800to1400_s,dy1400to2300_s,dy2300to3500_s,dy3500to4500_s,dy4500to6000_s,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]#inclmu15

raw_samples50m120 = [dy50to120,dy120to200,dy200to400,ttbar_pow,ww_incl,wz,zz_incl,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]
raw_samples120m = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,wz,zz_incl,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]

#raw_samples800m2300 = [dy800to1400_s,dy1400to2300,ttbar_pow]
#raw_samples400m2300 = [dy400to800_s,dy800to1400_s,dy1400to2300,ttbar_pow]

#raw_samples = [dy50to120, ttbar_pow]
#raw_samples = [dy50, ttbar, ww_incl, zz_incl, wz, dy50to120, ttbar_pow, dy50_startup, ttbar_startup]
#raw_samples = [zmumu, dy120_c1, dy200_c1, dy500_c1, dy1000_c1, ttbar, inclmu15]

# All MC samples
# lumi
#mcsum_lumi.prepare_histos_sum(raw_samples,lumiCD)
#mc_samples_sum_lumi = [nm1entry(sample,False,lumiBCD) for sample in raw_samples]
#for mc_sample in mc_samples_sum_lumi:
#    exec '%s = mc_sample' % mc_sample.name
# ref
#mcsum_ref.prepare_histos_sum(raw_samples, nolumi) 
#mc_samples_sum_ref = [nm1entry(sample,False,nolumi) for sample in raw_samples]
#for mc_sample in mc_samples_ref:
#    exec '%s = mc_sample' % mc_sample.name
#
# m > 50 GeV
# lumi
#mc50m_lumi.prepare_histos_sum(raw_samples50m, lumiBCD)
#mc_samples_50m_lumi = [nm1entry(sample,False,lumiBCD) for sample in raw_samples50m]
#for mc_sample in mc_samples_50m_lumi:
#    exec '%s = mc_sample' % mc_sample.name
# ref
#mc50m_ref.prepare_histos_sum(raw_samples50m, nolumi)
#mc_samples_50m_ref = [nm1entry(sample,False,nolumi) for sample in raw_samples50m]
#for mc_sample in mc_samples_ref:
#    exec '%s = mc_sample' % mc_sample.name
#
# 50 < m < 120 GeV
# lumi
#mc50m120_lumi.prepare_histos_sum(raw_samples50m120, lumiBCD)
mc_samples_50m120 = [nm1entry(sample,False,lumiBCD) for sample in raw_samples50m120]
for mc_sample in mc_samples_50m120:
    exec '%s = mc_sample' % mc_sample.name
# ref
#mc50m120_ref.prepare_histos_sum(raw_samples50m120, nolumi)
#mc_samples_ref = [nm1entry(sample,False,nolumi) for sample in raw_samples50m120]
#for mc_sample in mc_samples_ref:
#    exec '%s = mc_sample' % mc_sample.name
#
# M > 120 GeV
# lumi
#mc120m_lumi.prepare_histos_sum(raw_samples120m, lumiBCD)
mc_samples_120m = [nm1entry(sample,False,lumiBCD) for sample in raw_samples120m]
for mc_sample in mc_samples_120m:
    exec '%s = mc_sample' % mc_sample.name
# ref
#mc120m_ref.prepare_histos_sum(raw_samples120m, nolumi)
#mc_samples120m_ref = [nm1entry(sample,False,nolumi) for sample in raw_samples120m]
#for mc_sample in mc_samples120m_ref:
#    exec '%s = mc_sample' % mc_sample.name

## 800 < m < 2300 GeV
## lumi
#mc800m2300_lumi.prepare_histos_sum(raw_samples800m2300, lumiCD)
#mc_samples800m2300_lumi = [nm1entry(sample,False,lumiCD) for sample in raw_samples800m2300]
#for mc_sample in mc_samples800m2300_lumi:
#    exec '%s = mc_sample' % mc_sample.name
## ref
#mc800m2300_ref.prepare_histos_sum(raw_samples800m2300, nolumi)
#mc_samples800m2300_ref = [nm1entry(sample,False,nolumi) for sample in raw_samples800m2300]
#for mc_sample in mc_samples800m2300_ref:
#    exec '%s = mc_sample' % mc_sample.name
#
## 400 < m < 2300 GeV
## lumi
#mc400m2300_lumi.prepare_histos_sum(raw_samples400m2300, lumiCD)
#mc_samples800m2300_lumi = [nm1entry(sample,False,lumiCD) for sample in raw_samples800m2300]
#for mc_sample in mc_samples800m2300_lumi:
#    exec '%s = mc_sample' % mc_sample.name
## ref 
#mc400m2300_ref.prepare_histos_sum(raw_samples400m2300, nolumi)
#mc_samples800m2300_ref = [nm1entry(sample,False,nolumi) for sample in raw_samples800m2300]
#for mc_sample in mc_samples800m2300_ref:
#    exec '%s = mc_sample' % mc_sample.name
#
mass_ranges = [
#    ('50m',    (50, 1e9)),
#    ('70m',    (70, 1e9)),
#    ('70m110',  ( 70, 110)),
#    ('200m',    (200, 1e9)),
#    ('120m_BCD',   (120, 1e9)),
#    ('120m_CD',   (120, 1e9)),
#    ('60m120_BCD',  ( 60, 120)),
#    ('60m120_CD',  ( 60, 120)),
#    ('120m200', (120, 200)),
#    ('200m400', (200, 400)),
#    ('400m800', (400, 800)),
#    ('800m1400',(800, 1400)),
#    ('120m1400',(120,1400)),
#    ('all_lumi',  ( 120, 1e9)),
#    ('all_ref',  ( 120, 1e9)),
#    ('1400m2300',(1400, 2300)),
#    ('800m2300', (800,2300)),
#    ('400m2300', (400,2300)),
#    ('zpsi5000',(0,1e9)),
#    ('zpsi5000_m1TeV',(0,1e3)),
#    ('zpsi5000_1m3TeV',(1e3,3e3)),
#    ('zpsi5000_3mTeV',(3e3,1e9)),
     ('60m120',(60,120)),
#     ('60m',(60,2000)),
     ('120m',(120,2500)),
    ]

to_use = {
    
#    '60m120':  [data, dy50to120, dy50],#mcsum], #, zmumu, ttbar],
#    '70m110':  [data, dy50to120, ttbar_pow, mcsum],
#    '60m120':  [dataCD, mc50m120_lumi], #, zmumu, ttbar],
#    '120m200': [data, mcsum, dy50, ttbar, ww_incl, zz_incl, wz],#mcsum], #, dy120_c1, ttbar],
#    '200m400': [data, mcsum], #, dy200_c1, ttbar],
#    '400m600': [data, mcsum], #, dy200_c1, ttbar],
#    '200m':    [data, mcsum, dy50, ttbar, ww_incl, zz_incl, wz],#mcsum], #, dy500_c1, ttbar],
#    '50m':    [data, dy50to120],
#    '70m':    [data, dy50to120]
#    '120m_BCD':   [dataB,dataCD,mc120m_ref],
#    '120m':   [dataCD,mc120m_lumi],
#    '60m120':   [dataCD,mc50m120_lumi],
#    '60m120':   [dataBCD,mc_samples_50m120_lumi],
    '60m120':   [mc_samples_50m120,dataBCD],
#    '60m':      [dataCD,mc50m_lumi],
#    '120m':   [dataCD,mc120m_ref],
#    '120m':   [dataBCD,mc_samples_120m_lumi],
    '120m':   [mc_samples_120m,dataBCD],
#    '60m120':   [dataCD,mc50m120_ref],
#    '60m120_BCD':  [dataB,dataCD, mc50m120_ref], #, zmumu, ttbar],#powheg
#    'all_lumi': [dataCD,mc120m_lumi],
#    'all_ref': [dataCD,mc120m_ref],
#    '60m120_CD':  [dataCD, dy50to120,ttbar_pow],#mc50m120_lumi, zmumu, ttbar],#powheg
#    '120m200': [dataCD,DY120to200Powheg],
#    '200m400': [dataCD,DY200to400Powheg],
#    '400m800': [dataCD,DY400to800Powheg],
#    '800m1400': [dataCD,DY800to1400Powheg],
#    '1400m2300': [dataCD,dy1400to2300], # dataB,dataCD have no events and returns errors?
#    '800m2300': [dataCD,mc800m2300_lumi],
#    '400m2300': [dataCD,mc400m2300_lumi],
#    'zpsi5000':[zpsi5000],
#    'zpsi5000_m1TeV':[zpsi5000],
#    'zpsi5000_1m3TeV':[zpsi5000],
#    'zpsi5000_3mTeV':[zpsi5000],
    }

styles = {
#    'data':      (ROOT.kBlack,     -1),
    'dataB':      (ROOT.kBlack,     -1),
    'dataCD':      (ROOT.kBlack,     -1),
    'dataBCD':      (ROOT.kBlack,     -1),
#    'zmumu':     (ROOT.kRed,     3001),
#    'dy50':      (ROOT.kBlue,     1001),
#    'dy50_startup':      (ROOT.kViolet+2,     3001),
#    'dy200_c1':  (ROOT.kRed,     1001),
#    'dy500_c1':  (ROOT.kRed,     1001),
#    'dy1000_c1': (ROOT.kBlue,    1001),
#    'ttbar':     (ROOT.kGreen+2, 3005),
#    'ttbar_startup':     (ROOT.kGreen+2, 3005),
#    'ww_incl':    (ROOT.kBlue,    3004),
#    'zz_incl':    (62,            3006),
#    'wz' :       (64,            3007),
#    'inclmu15':  (28,            3002),
    'mc50m120_lumi':     (ROOT.kGreen+2, 1001),
    'mc50m120_ref':     (ROOT.kGreen+2, 1001),
    'mc50m_lumi':     (ROOT.kGreen+2, 1001),
    'mc50m_ref':     (ROOT.kGreen+2, 1001),
    'mc120m_lumi':     (ROOT.kGreen+2, 1001),
    'mc120m_ref':     (ROOT.kGreen+2, 1001),
    'mc800m2300_lumi':     (ROOT.kGreen+2, 1001),
    'mc800m2300_ref':     (ROOT.kGreen+2, 1001),
    'mc400m2300_lumi':     (ROOT.kGreen+2, 1001),
    'mc400m2300_ref':     (ROOT.kGreen+2, 1001),
    'dy50to120':      (ROOT.kGreen+2,     1001),
    'DY120to200Powheg':  (ROOT.kGreen+2, 1001),
    'DY200to400Powheg':  (ROOT.kGreen+2, 1001),
    'DY400to800Powheg':  (ROOT.kGreen+2, 1001),
    'DY800to1400Powheg': (ROOT.kGreen+2, 1001),
    'dy1400to2300':      (ROOT.kGreen+2, 1001),
    'ttbar_pow':     (ROOT.kGreen+2, 1001),
    'mcsum_lumi':    (ROOT.kGreen+2, 1001),
    'mcsum_ref':    (ROOT.kGreen+2, 1001),
    'zpsi5000':    (ROOT.kBlue, 1001),
    'zpsi5000_m1TeV':    (ROOT.kBlue, 1001),
    'zpsi5000_1m3TeV':    (ROOT.kBlue, 1001),
    'zpsi5000_3mTeV':    (ROOT.kBlue, 1001),
    'mc_samples_50m120': (ROOT.kGreen+2,1001),
    'mc_samples_120m': (ROOT.kGreen+2,1001),
    }

ymin = {
#    '70m110':  0.8,
#    '120m200': 0.6,
#    '200m400': 0.85,
#    '200m':    0.6,
#    '50m':    0.6,
#    '70m':    0.6,
#    '60m120_CD':  0.8,
    '60m':  0.6,
    '60m120':  0.7,
#    '60m120_BCD':  0.8,
#    '120m_BCD':    0.5,
#    '120m_CD':    0.5,
    '120m':    0.5,
    '120m200': 0.5,
    '200m400': 0.5,
    '400m800': 0.5,
    '800m1400': 0.5,
    'all_lumi':0.5,
    'all_ref':0.5,
    'zpsi5000':0.9,
    'zpsi5000_m1TeV':0.9,
    'zpsi5000_1m3TeV':0.9,
    'zpsi5000_3mTeV':0.9,
#    '120m1400':0.5,
#    '1400m2300':0.5,
#    '800m2300': 0.5,
#    '400m2300': 0.5,
    }
#global_ymin = 0.
global_ymin = None


def table_wald(entry, mass_range):
    print entry.name
    hnum = entry.histos['NoNo']
    mlo, mhi = mass_range
    num = get_integral(hnum, mlo, mhi, integral_only=True, include_last_bin=False)
    print 'mass range %5i-%5i:' % mass_range
    print 'numerator:', num
    print '%20s%20s%21s%25s%26s' % ('cut', 'denominator', 'efficiency', '- 68% CL-CP +','68% CL-Wald')
    for nminus1 in nminus1s:
        hden = entry.histos[nminus1]
        den = get_integral(hden, mlo, mhi, integral_only=True, include_last_bin=False)
        if num==0 and den==0:
            eff = 0
            errw = 0
        else:
            eff = float(num)/den
            errw = wald_binomial_err2(eff, 1, den)**0.5
        ecp,lcp,hcp = clopper_pearson(num, den)
        print '%20s%20f%20f%15f%15f%23f' % (nminus1, den, eff, eff-lcp, hcp-eff, errw)
        print '%20s%20f%20f%15f%15f%15f%16f' % (nminus1, den, eff, lcp, hcp, eff-errw, eff+errw)
        print ' '
    print '---------------------------------------------'

def print_wald_eff(y,eyl,eyh):
    print 'MC Sum'
    for i,nminus1 in enumerate(nminus1s):
        print nminus1, y[i], eyl[i], eyh[i]

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
        print '%20s%20f%20f%15f%15f' % (nminus1, den, e, e-l, h-e)

ROOT.gStyle.SetTitleX(0.25)
ROOT.gStyle.SetTitleY(0.50)

for name, mass_range in mass_ranges:
    pretty_name = pretty[name]
    print name, pretty_name

    lg = ROOT.TLegend(0.25, 0.21, 0.91, 0.44)
    lg.SetTextSize(0.03)
    lg.SetFillColor(0)
    lg.SetBorderSize(1)
    
    same = 'A'
    effs = []
    
    for entry in to_use[name]:

        l = len(nminus1s)
        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
            table(entry, mass_range)
            color, fill = styles[entry.name]
            nminus1_num = ROOT.TH1F('num', '', l, 0, l)
            nminus1_den = ROOT.TH1F('den', '', l, 0, l)
            hnum = entry.histos['NoNo']
            num = get_integral(hnum, *mass_range, integral_only=True, include_last_bin=False)
            for i,nminus1 in enumerate(nminus1s):
                hden = entry.histos[nminus1]
                den = get_integral(hden, *mass_range, integral_only=True, include_last_bin=False)
                nminus1_num.SetBinContent(i+1, num)
                nminus1_den.SetBinContent(i+1, den)
            eff,ycp,eylcp,eyhcp = binomial_divide(nminus1_num, nminus1_den)
        else:
            color = ROOT.kGreen+2
            fill = 1001
            nMC = len(entry)
            num = []
            pw = []
            den = []
            for i, mc in enumerate(entry):
                #table_wald(mc, mass_range)
                den.append([])
                hnum = mc.histos['NoNo']
                num.append(get_integral(hnum, *mass_range, integral_only=True, include_last_bin=False))
                pw.append(mc.partial_weight * lumiCD)
                for j, nminus1 in enumerate(nminus1s):
                    hden = mc.histos[nminus1]
                    den[i].append(get_integral(hden, *mass_range, integral_only=True, include_last_bin=False))
                    #if den[i][j]==0 and num[i]==0:
                    #    print mc.name, nminus1
            nm1Temp = ROOT.TH1F('temp','',l,0,l)
            eff,y,eyl,eyh = eff_wald(num, den, l, nMC, pw, nm1Temp)
            print_wald_eff(y,eyl,eyh)
        eff.SetTitle(pretty_name)
        eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin[name], 1.01)
        eff.GetXaxis().SetTitle('cut')
        eff.GetYaxis().SetLabelSize(0.027)
        eff.GetYaxis().SetTitle('n-1 efficiency')
        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
            draw = 'P'
            eff.SetLineColor(color)
            eff.SetMarkerStyle(20)
            eff.SetMarkerSize(1.05)
            eff.SetMarkerColor(color)
            #lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
            lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')
        else:
            draw = '2'
            eff.SetLineColor(color)
            eff.SetFillColor(color)
            eff.SetFillStyle(fill)
            #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            lg.AddEntry(eff, 'Simulation', 'LF')
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
