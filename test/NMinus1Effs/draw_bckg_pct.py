#!/usr/bin/env python

from pprint import pprint
import os,sys
from array import array
import numpy as np
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadRightMargin(0.03)
ROOT.gStyle.SetPadLeftMargin(0.12)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.TH1.AddDirectory(0)
outfile = ROOT.TFile("test.root","recreate")


psn = 'plots/paper2016'
ps = plot_saver(psn, size=(600,600), log=True, pdf=True, name='bckg_plot')

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
    'NoNo'
    ]

pretty = {
    'NoNo': 'All Selection',
    'NoPt': 'No p_{T} > 53 GeV',
    'NoDB': 'No |dxy| < 0.2',
    'NoIso': 'No rel. tk. iso.',
    'NoTkLayers': 'No # tk lay > 5',
    'NoPxHits': 'No # px hits > 0',
    'NoMuHits': 'No # mu hits > 0',
    'NoMuMatch': 'No # matched stations > 1',
    'NoVtxProb': 'No #chi^{2} #mu#mu vtx < 20',
    'NoB2B': 'No back-to-back',
    'NoDptPt': 'No dpT/pT',
    'NoTrgMtch': 'No HLT match',
    'ttbar_pow': 't#bar{t} powheg',
    'ww_incl': 'WW',
    'zz_incl': 'ZZ',
    'wz' : 'WZ',
    'dy50to120': 'DY#rightarrow#mu#mu 50 < m < 120 GeV',
    'dy120to200': 'DY#rightarrow#mu#mu 120 < m < 200 GeV',
    'dy200to400': 'DY#rightarrow#mu#mu 200 < m < 400 GeV',
    'dy400to800': 'DY#rightarrow#mu#mu 400 < m < 800 GeV',
    'dy800to1400': 'DY#rightarrow#mu#mu 800 < m < 1400 GeV',
    'dy1400to2300': 'DY#rightarrow#mu#mu 1400 < m < 2300 GeV',
    'dy2300to3500': 'DY#rightarrow#mu#mu 2300 < m < 3500 GeV',
    'dy3500to4500': 'DY#rightarrow#mu#mu 3500 < m < 4500 GeV',
    'dy4500to6000': 'DY#rightarrow#mu#mu 4500 < m < 6000 GeV',
    'tWtop': 'tW^{-}',
    'tWantitop': '#bar{t}W^{+}',
    'wjets': 'W + jets',
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
}

class nm1entry:
    def __init__(self, sample, is_data, lumi):
        if type(sample) == str:
            self.name = sample
            self.fn = self.make_fn(sample,is_data) if is_data else None
            self.lumi = lumi if is_data else None
        else:
            self.name = sample.name
            self.fn = self.make_fn(self.name,is_data)
            self.partial_weight = sample.partial_weight
        self.prepare_histos()
            
    def make_fn(self, name, is_data):
        if is_data==True:
            return 'data/ana_nminus1_%s.root' %name
        else:
            return 'mc/ana_nminus1_%s.root' % name
    
    def prepare_histos(self):
        self.histos = {}
        if self.fn is not None:
            f = ROOT.TFile(self.fn)
            for nminus1 in nminus1s:
                self.histos[nminus1] = f.Get(nminus1).Get('DimuonMassVertexConstrained').Clone()#DileptonMass

lumiBCD=2800.

from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
#raw_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,wjets,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]
raw_samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200]
use_samples = [qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets,tWtop,tWantitop,ww_incl,zz_incl,wz,ttbar_pow,dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000]

#mc_samples = [nm1entry(sample,False,lumiBCD) for sample in use_samples]
mc_samples = [nm1entry(sample,False,lumiBCD) for sample in reversed(raw_samples)]
for mc_sample in mc_samples:
    exec '%s = mc_sample' % mc_sample.name

to_use = {
    'NoPt':[mc_samples],
    'NoDB':[mc_samples],
    'NoIso':[mc_samples],
    'NoTkLayers':[mc_samples],
    'NoPxHits':[mc_samples],
    'NoMuHits':[mc_samples],
    'NoMuMatch':[mc_samples],
    'NoVtxProb':[mc_samples],
    'NoB2B':[mc_samples],
    'NoDptPt':[mc_samples],
    'NoTrgMtch':[mc_samples],
    'NoNo':[mc_samples],
}
styles = {
#    'sample':    (color, draw/fill style),
    #'dataB':      (ROOT.kBlack,     -1),
    #'dataCD':      (ROOT.kBlack,     -1),
    #'dataBCD':      (ROOT.kBlack,     -1),
    #'dy50to120':(ROOT.kGreen, 1001),
    #'dy120to200':(ROOT.kGreen-4, 1001),
    #'dy200to400':(ROOT.kGreen+1, 1001),
    #'dy400to800':(ROOT.kGreen-7, 1001),
    #'dy800to1400':(ROOT.kGreen-3, 1001),
    #'dy1400to2300':(ROOT.kGreen+2, 1001),
    #'dy2300to3500':(ROOT.kGreen-9, 1001),
    #'dy3500to4500':(ROOT.kGreen-6, 1001),
    #'dy4500to6000':(ROOT.kGreen-2, 1001),
    'dy50to120':(ROOT.kGreen+2, 1001),
    'dy120to200':(ROOT.kGreen+2, 1001),
    'dy200to400':(ROOT.kGreen+2, 1001),
    'dy400to800':(ROOT.kGreen+2, 1001),
    'dy800to1400':(ROOT.kGreen+2, 1001),
    'dy1400to2300':(ROOT.kGreen+2, 1001),
    'dy2300to3500':(ROOT.kGreen+2, 1001),
    'dy3500to4500':(ROOT.kGreen+2, 1001),
    'dy4500to6000':(ROOT.kGreen+2, 1001),
    'dy6000_s':(ROOT.kGreen+2, 1001),
    'ttbar_pow':(ROOT.kBlue,1001),
    'ww_incl':(ROOT.kOrange,1001),
    'zz_incl':(ROOT.kOrange,1001),
    'wz':(ROOT.kOrange,1001),
    'tWtop':(ROOT.kYellow,1001),
    'tWantitop':(ROOT.kYellow,1001),
    'wjets':(ROOT.kViolet,1001),
    'qcd50to80':(ROOT.kViolet,1001),
    'qcd80to120':(ROOT.kViolet,1001),
    'qcd120to170':(ROOT.kViolet,1001),
    'qcd170to300':(ROOT.kViolet,1001),
    'qcd300to470':(ROOT.kViolet,1001),
    'qcd470to600':(ROOT.kViolet,1001),
    'qcd600to800':(ROOT.kViolet,1001),
    'qcd800to1000':(ROOT.kViolet,1001),
    'qcd1000to1400':(ROOT.kViolet,1001),
    'qcd1400to1800':(ROOT.kViolet,1001),
    'qcd1800to2400':(ROOT.kViolet,1001),
    'qcd2400to3200':(ROOT.kViolet,1001),
    'qcd3200':(ROOT.kViolet,1001),
}

mass_range = [60,80,100,120]
maxX = 5000
bin_width = 50
nBins = maxX/bin_width
for i in xrange(3,nBins):
    ibin = i*bin_width
    mass_range.append(ibin)
nmc = len(mc_samples)
nnm1 = len(nminus1s)
nbins = len(mass_range)-1
sdata = (nnm1,nbins,nmc)
stotals = (nnm1,nbins)
data = np.zeros(sdata)
totals = np.zeros(stotals)
# - Make data
#   - data matrix 
#     3 dimensions (nminus1,mass_bin_lowEdge,mc_sample)
#   - totals matrix
#     2 dimensions (nminus1,mass_bin_lowEdge)
for nm1,nminus1 in enumerate(nminus1s):
    pretty_name = pretty[nminus1]
    print(nminus1,pretty_name)
    for entry in to_use[nminus1]: # single item (mc_samples)
        for ibin in range(len(mass_range)):
            if ibin == nbins: continue
            mlow = mass_range[ibin]
            mhigh = mass_range[ibin+1]
            total = 0
            for imc,mc in enumerate(entry):
                hmc = mc.histos[nminus1]
                count = get_integral(hmc,mlow,mhigh,integral_only=True,include_last_bin=False)
                count = count*sample.partial_weight
                data.itemset((nm1,ibin,imc),count)
                total = total+count
            totals.itemset((nm1,ibin),total)

# Make the Histograms
for nm1,nminus1 in enumerate(nminus1s):
    pretty_name = pretty[nminus1]
    lg = ROOT.TLegend(0.65,0.65,0.95,0.85)
    stack = ROOT.THStack('hs','')
    for entry in to_use[nminus1]:
        for mc,sample in enumerate(entry):
            mc_hist = ROOT.TH1F('mcHist','',nbins,array('f',mass_range))
            mc_hist.SetMaximum(1.)
            for ibin in range(len(mass_range)):
                if ibin == nbins: continue
                binContent = data[nm1,ibin,mc]/totals[nm1,ibin]
                mc_hist.SetBinContent(ibin,binContent)
            color,fill = styles[sample.name]
            mc_hist.SetLineColor(color)
            mc_hist.SetFillColor(color)
            stack.Add(mc_hist)
            if sample.name == 'dy50to120':
                lg.AddEntry(mc_hist,"Drell-Yan",'F' )
            elif sample.name == 'ttbar_pow':
                lg.AddEntry(mc_hist,"t#bar{t}","F")
            elif sample.name == 'ww_incl':
                lg.AddEntry(mc_hist,"DiBoson","F")
            elif sample.name == 'qcd80to120':
                lg.AddEntry(mc_hist,"QCD & W+jets","F")
            elif sample.name == 'tWtop':
                lg.AddEntry(mc_hist,"Single Top","F")
    stack.SetMaximum(1.)
    stack.SetMinimum(1.E-4)
    stack.Draw('hist')
    stack.SetTitle(pretty_name+' Background Percentage')
    stack.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    stack.GetYaxis().SetTitle('percentage/%i GeV'%bin_width)
    stack.GetYaxis().SetTitleOffset(1.4)
    lg.Draw()
    ps.save(nminus1+'_bckg_pct',pdf=True,log=True)
