#!/usr/bin/env python

# Imports and style settings to be applied to all scripts
from pprint import pprint
import sys, os
from array import array
import numpy as np
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()



#tt = 'tt pow'
tt = 'tt lep'

# lumi = 20409.950
# period = 'BF'
# period_1 = 'RunB - RunF'

# lumi = 16786.820
# period = 'GH' 
# period_1 = 'RunG - RunH'

lumi = 36295.39
period = '' 

folder = 'LeptonPt'
axisX = 'p_{T}(#mu) [GeV]'
# folder = 'DileptonPt'
# folder = 'NTkLayers'
# folder = 'RelIsoSumPt'
# folder = 'NPxHits'
# folder = 'DimuonMuonPtErrOverPt'
# folder = 'DimuonMassVtx_chi2'
# axisX = 'Vtx #chi^{2}'

# folder = 'DileptonMass'
# folder = 'DileptonMass_bb'
# folder = 'DileptonMass_be'
# folder = 'DimuonMassVertexConstrained'
# axisX = 'm(#mu#mu) [GeV]'
# folder = 'DimuonMassVertexConstrained_bb'
# folder = 'DimuonMassVertexConstrained_be'

styles = {
#    'sample':    (color, draw/fill style),
    'data':      (ROOT.kBlack,     -1),
    'dataCD':      (ROOT.kBlack,     -1),
    'dataBCD':      (ROOT.kBlack,     -1),


#                                 return 'DY #rightarrow #tau#tau', 5
#                         else:
#                                 return  'DY #rightarrow #mu#mu', 3
 
  
     
          
    'dyInclusive50': (5,1001),
    'dy50to120':(3, 1001),
    'dy120to200':(3, 1001),
    'dy200to400':(3, 1001),
    'dy400to800':(3, 1001),
    'dy800to1400':(3, 1001),
    'dy1400to2300':(3, 1001),
    'dy2300to3500':(3, 1001),
    'dy3500to4500':(3, 1001),
    'dy4500to6000':(3, 1001),
#    'dy6000':(ROOT.kGreen+2, 1001),
    'ttbar_lep':(4,1001),
    'ttbar_lep50to500':(4,1001),
    'ttbar_lep50to800':(4,1001),
    'ttbar_lep_500to800':(4,1001),
    'ttbar_lep_800to1200':(4,1001),
    'ttbar_lep_1200to1800':(4,1001),
    'ttbar_lep_1800':(4,1001),
    'ttbar_pow':(ROOT.kBlue,1001),
    'WWinclusive':(ROOT.kRed-10,1001),
    'WW200to600':(ROOT.kRed-10,1001),
    'WW600to1200':(ROOT.kRed-10,1001),
    'WW1200to2500':(ROOT.kRed-10,1001),
    'WW2500':(ROOT.kRed-10,1001),
    'ZZ':(2,1001),
    'ZZ_ext':(2,1001),
    'WZ':(ROOT.kRed+3,1001),
    'WZ_ext':(ROOT.kRed+3,1001),
    'tW':(ROOT.kBlue+2,1001),
    'Wantitop':(ROOT.kBlue+2,1001),
    'Wjets':(7,1001),
    'qcd80to120':(20,1001),
    'qcd120to170':(20,1001),
    'qcd170to300':(20,1001),
    'qcd300to470':(20,1001),
    'qcd470to600':(20,1001),
    'qcd600to800':(20,1001),
    'qcd800to1000':(20,1001),
    'qcd1000to1400':(20,1001),
    'qcd1400to1800':(20,1001),
    'qcd1800to2400':(20,1001),
    'qcd2400to3200':(20,1001),
    'qcd3200':(20,1001),


}

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

tightnm1 = [
    #'TiPt',
    'TiDB',
    'TiGlbChi2',
    'TiIso',
    'TiTkLayers',
    'TiPxHits',
    'TiMuHits',
    'TiMuMatch',
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
    'NoNo': 'All Selection',
    #'data': 'Data, %.1f fb^{-1}',
    'data': 'Data, %.1f fb^{-1}, MuonOnly',
    'dataB': 'Data RunB, %.1f fb^{-1}, MuonOnly',
    'dataCD': 'Data RunC+D, %.1f fb^{-1}, MuonOnly',
    'dataBCD': 'Data, %.1f fb^{-1}, 2015 MuonOnly',
    'all': 'm > 60 GeV',
    '120m': 'm > 120 GeV',
    '500m': 'm > 500 GeV',
    '60m120': '60 < m < 120 GeV',
    'zpsi5000': 'Z\'_{#psi}, M=5000 GeV',
    'dy50to120': 'DY#rightarrow#mu#mu 50 < m < 120 GeV',
    'dy120to200': 'DY#rightarrow#mu#mu 120 < m < 200 GeV',
    'dy200to400': 'DY#rightarrow#mu#mu 200 < m < 400 GeV',
    'dy400to800': 'DY#rightarrow#mu#mu 400 < m < 800 GeV',
    'dy800to1400': 'DY#rightarrow#mu#mu 800 < m < 1400 GeV',
    'dy1400to2300': 'DY#rightarrow#mu#mu 1400 < m < 2300 GeV',
    'dy2300to3500': 'DY#rightarrow#mu#mu 2300 < m < 3500 GeV',
    'dy3500to4500': 'DY#rightarrow#mu#mu 3500 < m < 4500 GeV',
    'dy4500to6000': 'DY#rightarrow#mu#mu 4500 < m < 6000 GeV',
    'dy6000': 'DY#rightarrow#mu#mu m > 6000 GeV',
    'ttbar_lep': 't#bar{t} lep',
    'ttbar_lep50to500': 't#bar{t} lep 50 < m < 500 GeV',
    'ttbar_lep50to800': 't#bar{t} lep 50 < m < 500 GeV',
    'ttbar_lep_800to1200': 't#bar{t} lep',
    'ttbar_lep_1200to1800': 't#bar{t} lep',
    'ttbar_lep1800toInf': 't#bar{t} lep',
    'ttbar_pow': 't#bar{t} powheg',
    'ww_incl': 'WW',
    'zz_incl': 'ZZ',
    'wz' : 'WZ',
    'tWtop': 'tW^{-}',
    'tWantitop': '#bar{t}W^{+}',
    'wjets': 'W + jets',
    'qcd50to80': 'QCD 50 < p_{T}(#mu) < 80 GeV',
    'qcd80to120': 'QCD 80 < p_{T}(#mu) < 120 GeV',
    'qcd120to170': 'QCD 120 < p_{T}(#mu) < 170 GeV',
    'qcd170to300': 'QCD 170 < p_{T}(#mu) < 300 GeV',
    'qcd300to470': 'QCD 300 < p_{T}(#mu) < 470 GeV',
    'qcd470to600': 'QCD 470 < p_{T}(#mu) < 600 GeV',
    'qcd600to800': 'QCD 600 < p_{T}(#mu) < 800 GeV',
    'qcd800to1000': 'QCD 800 < p_{T}(#mu) < 1000 GeV',
    'qcd1000to1400': 'QCD 1000 < p_{T}(#mu) < 1400 GeV',
    'qcd1400to1800': 'QCD 1400 < p_{T}(#mu) < 1800 GeV',
    'qcd1800to2400': 'QCD 1800 < p_{T}(#mu) < 2400 GeV',
    'qcd2400to3200': 'QCD 2400 < p_{T}(#mu) < 3200 GeV',
    'qcd3200': 'QCD p_{T}(#mu) < 3200 GeV',
    'inclmu15': 'QCD',
}

def wald_binomial_err2(p_hat, pw, den):
    return float(p_hat) * (1 - float(p_hat)) / (float(den)/pw) # den is pre-weighted number of MC events

def eff_wald(num, den, l, nMC, pw,temp):
    '''
    - Only used in NMinus1Effs/draw.py
    - Inputs
        - num = Numerator MC histogram (per NM1 eff)
        - den = Denominator MC histogram (per NM1 eff)
        - l = len(nminus1s)
        - nMC = len(entry), entry = [dy50to120, ...]
        - pw = [dy50to120.partial_weight*lumi, ...]
        - temp = ROOT.TH1F('temp','',l,0,l)
    - Output
        - MC N-1 Histogram with Wald-Approx. uncertainties
        - Array of effs and +/- uncertainties for each NM1 bin
    '''
    nbins = temp.GetNbinsX()
    xax = temp.GetXaxis()
    x = []
    y = []
    exl = []
    exh = []
    eyl = []
    eyh = []
    numTot = sum(e*f for e,f in zip(pw,num))
    for i in range(l): # i = NoX
        ibin = i+1
        x.append(xax.GetBinCenter(ibin))
        exl.append(xax.GetBinWidth(ibin)/2)
        exh.append(xax.GetBinWidth(ibin)/2)
        denTot_i = 0
        sum_err2_i = 0
        err2_i = []
        denTot_i = sum([b*c for b,c in zip(pw,[a[i] for a in den])])
        if numTot>denTot_i:
            print "in eff_wald(): numTot>denTot, setting eff=1 ",numTot,denTot_i
            y.append(1.)
        elif denTot_i!=0:
            y.append(numTot/denTot_i)
        else: 
            print "in eff_wald(): denTot_i==0 ?, setting eff=0",denTot_i
            y.append(0.)
        for mc in range(nMC): 
            #print float(num[mc]), den[mc][i]
            p_hat_i_mc = 0
            if num[mc] != 0:
                p_hat_i_mc = float(num[mc])/den[mc][i] # mc = mc sample, i = NoX
                err2_i.append((pw[mc]*den[mc][i]/denTot_i)**2 * wald_binomial_err2(p_hat_i_mc,1., den[mc][i]))
            else:
                p_hat_i_mc = 0
                err2_i.append(0)
            #print p_hat_i_mc
        sum_err2_i = sum(err2_i)
        eyh.append(sum_err2_i**0.5)
        eyl.append(sum_err2_i**0.5)
    eff = ROOT.TGraphAsymmErrors(len(x), *[array('d', obj) for obj in (x,y,exl,exh,eyl,eyh)])
    return eff, y, eyl, eyh


def print_wald_eff(y,eyl,eyh):
    print 'MC Sum'
    for i,nminus1 in enumerate(nminus1s):
        print nminus1, y[i], eyl[i], eyh[i]

def table_wald(entry,nminus1, mass_range):
    print entry.name, nminus1
    hnum = entry.histos['NoNo']
    hden = entry.histos[nminus1]
    print '%20s%27s%23s%20s%16s%25s%26s' % ('cut', 'mass range','numerator', 'denominator', 'efficiency', '- 68% CL-CP +','68% CL-Wald')
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
        	num = num*entry.partial_weight*lumi
        	den = den*entry.partial_weight*lumi
#         print '%20s%15i%15i%20f%20f%15f%15f%15f%23f'     % (nminus1, mlow, mhigh, num, den, eff, eff-lcp, hcp-eff,        errw)
#         print '%20s%15i%15i%20f%20f%15f%15f%15f%15f%16f' % (nminus1, mlow, mhigh, num, den, eff, lcp,     hcp,     eff-errw, eff+errw)
        print '%20s%15i%15i%20f%20f%15f%15f%15f%23f'     % (nminus1, mlow, mhigh, num, den, eff, eff-lcp, hcp-eff,        errw)
        print '%20s%15i%15i%20f%20f%15f%15f%15f%15f%16f' % (nminus1, mlow, mhigh, num, den, eff, lcp,     hcp,     eff-errw, eff+errw)

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
        print '%20s%15i%15i%20f%20f%20f%15f%15f' % (nminus1, mlow, mhigh, num, den, e, l, h)

# nm1entry
class nm1entry:
    def __init__(self, sample, is_data, lumi):
        if type(sample) == str:
            self.name = sample      
            self.fn = '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/NMinus1Effs/data/Run2016MuonsOnly%s/ana_nminus1_data.root' %period
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
        #return 'mc/ana_nminus1_%s.root' % name
        return '/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/NMinus1Effs/mc/ana_nminus1_%s.root' %  name       
        #return 'mc_TriggerScale/ana_nminus1_%s.root' % name
        
        
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s + ['NoNo']:
                self.histos[nminus1] = f.Get(nminus1).Get(folder).Clone()
#                self.histos[nminus1] = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()

    # This function isn't used anymore, but keep for posterity? cjsbad
    def prepare_histos_sum(self, samples, lumi):
        self.histos = {}
        for nminus1 in nminus1s + ['NoNo']:
            hs = []
            print '%20s%20s%21s%20s%20s' % ('cut', 'sampe name', 'partial weight', 'scale(ref)','scale(lumi)')
            for sample in samples:
                f = ROOT.TFile(self.make_fn(sample.name))
                #h = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()
                h = f.Get(nminus1).Get(folder).Clone()
                if lumi>0:
                    # scale to luminosity for comparision of single dataset to MC
                    h.Scale(sample.partial_weight * lumi) 
                if lumi<0:
                    # scale to reference cross section/Nevents for comparision of multiple datasets to MC
                    h.Scale(sample.partial_weight * refN / refXS)  
                hs.append(h)
            hsum = hs[0].Clone()
            for h in hs[1:]:
                hsum.Add(h)
            self.histos[nminus1] = hsum

