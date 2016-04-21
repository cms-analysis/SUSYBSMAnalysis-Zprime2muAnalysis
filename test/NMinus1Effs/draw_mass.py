#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

from pprint import pprint
import sys, os
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)
ROOT.TH1.AddDirectory(0)

outfile = ROOT.TFile("whargl.root","recreate")
iarp=0
mcN=0
drawDY = True
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
            self.fn = 'data/ana_nminus1_%s.root'%sample
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
        return 'mc/ana_nminus1_%s.root' % name
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                self.histos[nminus1] = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()#DileptonMass
                '''
                if nminus1=='NoVtxProb':
                    self.histos[nminus1] = f.Get(nminus1).Get('DileptonMass').Clone()
                else:
                    self.histos[nminus1] = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()#DileptonMass
                '''

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
#DYmc = nm1entry('DYmc',False,lumiBCD)
#nonDYmc = nm1entry('nonDYmc',False,lumiBCD)
#wjets = nm1entry('wjets',False,lumiBCD)

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
#raw_samples = [dy50to120,DY120to200Powheg,DY200to400Powheg,DY400to800Powheg,DY800to1400Powheg,dy1400to2300,dy2300to3500,DY3500to4500Powheg,dy4500to6000,ttbar_pow_s,ww_incl_s,zz_incl,wz,tWtop,tWantitop,wjets]#inclmu15,
#,qcd600to800,qcd120to170
#raw_samples = [dy50to120_s,dy120to200_s,dy200to400_s,dy400to800_s,dy800to1400_s,dy1400to2300_s,dy2300to3500_s,dy3500to4500_s,dy4500to6000_s,dy6000_s,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,inclmu15,zpsi5000,qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]

raw_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,zpsi5000,wjets]#dy6000_s,
use_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]#dy6000_s,

#dy_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000]#dy6000,

#nonDY_samples = [ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop,qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]#dy6000_s,

#EWKnoDY_samples = [ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop]#dy6000_s,

#noQCD_samples = [dy50to120_s,dy120to200_s,dy200to400_s,dy400to800_s,dy800to1400_s,dy1400to2300_s,dy2300to3500_s,dy3500to4500_s,dy4500to6000_s,ttbar_pow,ww_incl,zz_incl,wz,wjets,tWtop,tWantitop]

#QCD_samples = [qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200]

#diboson_samples = [ww_incl, zz_incl, wz]

#stop_samples = [tWtop, tWantitop]

#ttbar_samples = [ttbar_pow]

#wjets_samples = [wjets]

#noQCDwjets_samples = [dy50to120_s,dy120to200_s,dy200to400_s,dy400to800_s,dy800to1400_s,dy1400to2300_s,dy2300to3500_s,dy3500to4500_s,dy4500to6000_s,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop]

#QCDwjets_samples = [qcd50to80,qcd80to120,qcd170to300,qcd300to470,qcd470to600,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd3200,wjets]

#zpsi5000_samples = [zpsi5000]

#EWKnoDYwjets_samples = [ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop]

#refXS = dy50to120_s.cross_section
#refN = dy50to120_s.nevents

# All MC samples
# lumi
#mcsum_lumi.prepare_histos_sum(use_samples,lumiBCD)
mc_samples = [nm1entry(sample,False,lumiBCD) for sample in use_samples]
for mc_sample in mc_samples:
    exec '%s = mc_sample' % mc_sample.name

#dy = [nm1entry(sample,False,lumiBCD) for sample in dy_samples]
#for dy_sample in dy_samples:
#    exec '%s = dy_sample' % dy_sample.name
#noDY = [nm1entry(sample,False,lumiBCD) for sample in nonDY_samples]
#for nonDY_sample in nonDY_samples:
#    exec '%s = nonDY_sample' % nonDY_sample.name
#qcd = [nm1entry(sample,False,lumiBCD) for sample in QCD_samples]
#for QCD_sample in QCD_samples:
#    exec '%s = QCD_sample' % QCD_sample.name
#diboson = [nm1entry(sample,False,lumiBCD) for sample in diboson_samples]
#stop = [nm1entry(sample,False,lumiBCD) for sample in stop_samples]
#ttbar = [nm1entry(sample,False,lumiBCD) for sample in ttbar_samples]
#noQCD = [nm1entry(sample,False,lumiBCD) for sample in noQCD_samples]
#wjets = [nm1entry(sample,False,lumiBCD) for sample in wjets_samples]
#noQCDwjets = [nm1entry(sample,False,lumiBCD) for sample in noQCDwjets_samples]
#EWKnoDY = [nm1entry(sample,False,lumiBCD) for sample in EWKnoDY_samples]
#QCDwjets = [nm1entry(sample,False,lumiBCD) for sample in QCDwjets_samples]
#Zpsi5000 = [nm1entry(sample,False,lumiBCD) for sample in zpsi5000_samples]
#EWKnoDYwjets = [nm1entry(sample,False,lumiBCD) for sample in EWKnoDYwjets_samples]
#

# ref
#mcsum_ref.prepare_histos_sum(use_samples, nolumi) 
#mc_samples_ref = [nm1entry(sample,False,nolumi) for sample in use_samples]
#for mc_sample in mc_samples_ref:
#    exec '%s = mc_sample' % mc_sample.name

#DYmc.prepare_histos_sum(DYmc_list,lumiBCD)
#nonDYmc.prepare_histos_sum(nonDYmc_list,lumiBCD)

#bin_width = 20
#maxX = 2500
#minX = 60
#nBins = (maxX-minX)/bin_width
#mass_range = []
#for i in range(3,nBins):
#    ibin = i*bin_width
#    mass_range.append(ibin)
#print mass_range
mass_range = [60,120,180,240,320,500,1000,2500]

to_use = {
#   'sample':[MC,Data],
    #'NoPt':[DYmc,nonDYmc,dataCD],
    #'NoDB':[DYmc,nonDYmc,dataCD],
    #'NoIso':[DYmc,nonDYmc,dataCD],
    #'NoTkLayers':[DYmc,nonDYmc,dataCD],
    #'NoPxHits':[DYmc,nonDYmc,dataCD],
    #'NoMuHits':[DYmc,nonDYmc,dataCD],
    #'NoMuMatch':[DYmc,nonDYmc,dataCD],
    #'NoVtxProb':[DYmc,nonDYmc,dataCD],
    #'NoB2B':[DYmc,nonDYmc,dataCD],
    #'NoDptPt':[DYmc,nonDYmc,dataCD],
    #'NoTrgMtch':[DYmc,nonDYmc,dataCD],

    'NoPt':[mc_samples,dataBCD],
    'NoDB':[mc_samples,dataBCD],
    'NoIso':[mc_samples,dataBCD],
    'NoTkLayers':[mc_samples,dataBCD],
    'NoPxHits':[mc_samples,dataBCD],
    'NoMuHits':[mc_samples,dataBCD],
    'NoMuMatch':[mc_samples,dataBCD],
    'NoVtxProb':[mc_samples,dataBCD],
    'NoB2B':[mc_samples,dataBCD],
    'NoDptPt':[mc_samples,dataBCD],
    'NoTrgMtch':[mc_samples,dataBCD],

    #'NoPt':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoDB':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoIso':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoTkLayers':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoPxHits':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoMuHits':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoMuMatch':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoVtxProb':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoB2B':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoDptPt':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    #'NoTrgMtch':[EWKnoDY,Zpsi5000,noDY,dy,dataBCD],
    }

styles = {
#    'sample':    (color, draw/fill style),
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
    'DYmc': (ROOT.kGreen, 3001),
    'nonDYmc': (ROOT.kRed, 1001),
    'dy':(ROOT.kGreen,1001),
    }

yrange = {
#   'sample':    (ymin,ymax),
    'NoPt':      (0.00,1.01),
    'NoDB':      (0.95,1.001),
    'NoIso':     (0.60,1.01),
    'NoTkLayers':(0.95,1.001),
    'NoPxHits':  (0.80,1.001),
    'NoMuHits':  (0.80,1.001),
    'NoMuMatch': (0.80,1.005),
    'NoVtxProb': (0.90,1.001),
    'NoB2B':     (0.80,1.001),
    'NoDptPt':   (0.95,1.001),
    'NoTrgMtch': (0.90,1.001),
    }
#global_ymin = 0.
global_ymin = None

def table_wald(entry,nminus1, mass_range):
    print entry.name, nminus1
    hnum = entry.histos['NoNo']
    hden = entry.histos[nminus1]
    #print '%20s%27s%23s%20s%16s%25s%26s' % ('cut', 'mass range','numerator', 'denominator', 'efficiency', '- 68% CL-CP +','68% CL-Wald')
    for mbin in range(len(mass_range)):
        if mbin == (len(mass_range)-1): break
        mlow = mass_range[mbin]
        mhigh = mass_range[mbin+1] 
        num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        pcp,lcp,hcp = clopper_pearson(num, den)
        if num==0 and den==0:
            eff = 0
            errw = 0
        else:
            eff = num/den
            if (eff*(1-eff))<0:
                print "what is this"
                print 'NM1, sample, mlow, mhigh, num, den'
                print nminus1, entry.name, mlow, mhigh, num, den
                eff = 1
                errw = 0
            else:
                errw = (eff*(1-eff)/den)**0.5
        if 'data' not in entry.name:
            num = num*entry.partial_weight*lumiBCD
            den = den*entry.partial_weight*lumiBCD
        #print '%20s%15i%15i%20f%20f%15f%15f%15f%23f'     % (nminus1, mlow, mhigh, num, den, eff, eff-lcp, hcp-eff,        errw)
        #print '%20s%15i%15i%20f%20f%15f%15f%15f%15f%16f' % (nminus1, mlow, mhigh, num, den, eff, lcp,     hcp,     eff-errw, eff+errw)
        print ' '
    print '---------------------------------------------'

def table(entry,nminus1, mass_range):
    print entry.name
    hnum = entry.histos['NoNo']
    hden = entry.histos[nminus1]
    print '%20s%27s%23s%20s%20s%22s' % ('cut', 'mass range', 'numerator', 'denominator', 'efficiency', '68% CL')
    for mbin in range(len(mass_range)):
        if mbin == (len(mass_range)-1): break
        mlow = mass_range[mbin]
        mhigh = mass_range[mbin+1] 
        num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
        den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
        e,l,h = clopper_pearson(num, den)
        print '%20s%15i%15i%20f%20f%20f%15f%15f' % (nminus1, mlow, mhigh, num, den, p_hat, p_hat_e, p_hat_e)

ROOT.gStyle.SetTitleX(0.25)
ROOT.gStyle.SetTitleY(0.50)

for nminus1 in nminus1s:
    pretty_name = pretty[nminus1]
    print nminus1, pretty_name
    lg = ROOT.TLegend(0.25, 0.21, 0.91, 0.44)
    lg.SetTextSize(0.03)
    lg.SetFillColor(0)
    lg.SetBorderSize(1)
    
    same = 'A'
    effs = []



    for entry in to_use[nminus1]: #,mass_range 

        #table(entry,nminus1, mass_range)

        l = len(mass_range)-1
        nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
        nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))

        if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
            #table_wald(entry,nminus1,mass_range)
            color, fill = styles[entry.name]
            hnum = entry.histos['NoNo']
            hden = entry.histos[nminus1]
            for mbin in range(len(mass_range)):
                if mbin == (len(mass_range)-1): continue
                mlow = mass_range[mbin]
                mhigh = mass_range[mbin+1]
                num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                nminus1_num.SetBinContent(mbin+1, num)
                nminus1_den.SetBinContent(mbin+1, den)
            eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
        else:
            #for a,mc in enumerate(entry):
            #    table_wald(mc,nminus1,mass_range)
            p_hats = []
            errsW = []
            x = []
            ex = []
            for mbin in range(len(mass_range)):
                if mbin == (len(mass_range)-1): continue
                numTot = 0
                denTot = 0
                err2sum = 0
                numsW = []
                densW = []
                err2s = []
                for i,mc in enumerate(entry):
                    hnum = mc.histos['NoNo']
                    hden = mc.histos[nminus1]
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False)
                    denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False)
                    if numInt<=denInt and denInt>0:
                        p_hat_mc = float(numInt)/denInt
                        err2 = p_hat_mc*(1-p_hat_mc)/denInt
                    elif numInt>denInt:
                        print "numInt>denInt"
                        print "NM1, sample, mlow, mhigh, num, den"
                        print nminus1, mc.name, mlow, mhigh, numInt, denInt
                        print
                        p_hat_mc = 1
                        err2 = 0
                    else:
                        if numInt!=0 or denInt!=0:
                            print "ELSE"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc.name, mlow, mhigh, numInt, denInt
                            print
                        p_hat_mc = 0
                        err2 = 0
                    numTot = numTot + numInt*mc.partial_weight
                    denTot = denTot + denInt*mc.partial_weight
                    numsW.append(numInt*mc.partial_weight)
                    densW.append(denInt*mc.partial_weight)
                    err2s.append(err2)
                p_hat = float(numTot)/denTot
                err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW,err2s))
                if err2sum<0:
                    err2sum = 0
                elif p_hat==0:
                    err2sum = 0
                p_hats.append(p_hat)
                errsW.append(err2sum**0.5)
                x.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                ex.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
            eff = ROOT.TGraphAsymmErrors(len(x), *[array('d',obj) for obj in (x,p_hats,ex,ex,errsW,errsW)])
        eff.SetTitle(pretty_name)
        ymin, ymax = yrange[nminus1]
        eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
        eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
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
            if len(entry)==9:
                draw = '2'
                eff.SetLineColor(ROOT.kGreen+2)
                eff.SetFillColor(ROOT.kGreen+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'DY MC 76X','LF')
            elif len(entry)==15:
                draw = '2'
                eff.SetLineColor(ROOT.kOrange+2)
                eff.SetFillColor(ROOT.kOrange+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK MC','LF')
            elif len(entry)==1:
                draw = '2'
                eff.SetLineColor(ROOT.kYellow+1)
                eff.SetFillColor(ROOT.kYellow+1)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'Z\'_{#psi} (5 TeV)','LF')
            elif len(entry)==10:
                draw = '2'
                eff.SetLineColor(ROOT.kViolet)
                eff.SetFillColor(ROOT.kViolet)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'QCD MC','LF')
            elif len(entry)==11:
                draw = '2'
                eff.SetLineColor(ROOT.kBlue+2)
                eff.SetFillColor(ROOT.kBlue+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'QCD + W+jets MC','LF')
            elif len(entry)==7:
                draw = '2'
                eff.SetLineColor(ROOT.kRed+2)
                eff.SetFillColor(ROOT.kRed+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK no DY MC','LF')
            elif len(entry)==6:
                draw = '2'
                eff.SetLineColor(ROOT.kRed+2)
                eff.SetFillColor(ROOT.kRed+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK','LF')
            elif len(entry)==len(use_samples):
                draw = '2'
                eff.SetLineColor(ROOT.kGreen+2)
                eff.SetFillColor(ROOT.kGreen+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'All Simulation','LF')
            else:
                draw = '2'
                eff.SetLineColor(ROOT.kBlue+2)
                eff.SetFillColor(ROOT.kBlue+2)
                eff.SetFillStyle(1001)
                lg.AddEntry(eff,'EWK+QCD MC','LF')
        #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
        draw += same
        eff.Draw(draw)
        effs.append(eff)
        same = ' same'
        outfile.cd()
        eff.Write("arp%d"%iarp)
        iarp+=1
    # end for entry in to_use[name]: # entry is a specific sample
    lg.Draw()
    ps.save(nminus1+'_mass')
    print
# end for name, mass_range in mass_bins:
