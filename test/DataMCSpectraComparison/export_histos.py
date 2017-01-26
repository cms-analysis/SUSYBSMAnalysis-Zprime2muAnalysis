#!/usr/bin/env python

# Script to export our histograms into the same format that HEEP uses
# so they can be drawn with a common script (plot_for_paper.C).

import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *

ROOT.TH1.AddDirectory(False)

mumu_mc_fn_base = 'mc/ana_datamc_%s.root'
mumu_data_fn = 'data/Run2015MuonsOnly/ana_datamc_data.root'
mumu_rebin_factor = 1
mumu_histogram = 'DimuonMassVertexConstrained'
mumu_scale = 48
#mumu_scale = 1.
emu_mc_fn_base = 'mc/ana_datamc_%s.root'
emu_data_fn = 'data/Run2015MuonsOnly/ana_datamc_data.root' #it was Run2015 only
emu_rebin_factor = 20
emu_histogram = 'DileptonMass'
emu_scale = 4641
add_heep = False
heep_fn = 'heep_massHists2011Mar13VarBinsDoubleRes.root'
heep_rebin_factor = None

mumu_rebin_roundto = None
# resolution = lambda x: 0.009332*x + 5.71e-5*x*x - 1.171e-9*x*x*x
#mumu_rebin_variable = [i for i in xrange(60, 255, 5)] + [260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 365, 380, 395, 410, 425, 440, 460, 480, 500, 520, 545, 570, 595, 625, 655, 690, 725, 765, 810, 855, 905, 960, 1025, 1095, 1175, 1265, 1370, 1490, 1630, 1795, 1990, 2230, 2500]
mumu_rebin_variable = [i for i in xrange(60, 150, 5)] + [i for i in xrange(150, 200, 10)] + [i for i in xrange(200, 600, 20)] + [i for i in xrange(600, 900, 30)] + [i for i in xrange(900, 1250, 40)]
# resolution 2x above
#mumu_rebin_variable = [i for i in xrange(60, 145, 5)] + [155, 165, 175, 185, 195, 205, 215, 225, 235, 250, 265, 280, 295, 315, 335, 355, 380, 405, 435, 465, 500, 540, 585, 635, 695, 765, 850, 950, 1070, 1220, 1410, 1660, 1995, 2500]
#mumu_rebin_variable = None
heep_rebin_roundto = 1

from array import array
def variable_rebin_and_normalize(h, bins, roundto):
    if bins is not None:
        bins = array('d', bins)
        h = h.Rebin(len(bins)-1, h.GetName() + '_rebin', bins)
    if roundto is not None:
        for i in xrange(1, h.GetNbinsX()+1):
            width = h.GetBinWidth(i)
            if width > roundto:
                h.SetBinContent(i, h.GetBinContent(i) * roundto / width)
    return h
    
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

# First mu+mu-, differential and cumulative.
for cumulative in (False, True):
    histos = {}
    for sample in samples:
        f = ROOT.TFile(mumu_mc_fn_base % sample.name)
        h = f.Histosstartup.Get(mumu_histogram).Clone()
        if mumu_rebin_variable and not cumulative:
            h = variable_rebin_and_normalize(h, mumu_rebin_variable, mumu_rebin_roundto)
        else:
            h.Rebin(mumu_rebin_factor)
        h.Scale(sample.partial_weight * mumu_scale)
        if cumulative:
            h = cumulative_histogram(h)
            h.SetName(h.GetName().replace('_cumulative_ge', ''))
        histos[sample.name] = h

    f = ROOT.TFile(mumu_data_fn)
    dataHist = f.Histosstartup.Get(mumu_histogram).Clone('data')
    if mumu_rebin_variable and not cumulative:
        dataHist = variable_rebin_and_normalize(dataHist, mumu_rebin_variable, mumu_rebin_roundto)
        dataHist.SetName('data')
    else:
        dataHist.Rebin(mumu_rebin_factor)
        
    if cumulative:
        dataHist = cumulative_histogram(dataHist)
        dataHist.SetName('data')
    
#    jetsHist = histos['inclmu15'].Clone('jetsHist')
    jetsHist = histos['wjets'].Clone('jetBkgHist')
#    jetsHist.Add(histos['wjets'])

    promptHist = histos['ttbar_pow'].Clone('prompt')
    promptHist.Add(histos['tWtop'])
    promptHist.Add(histos['tWantitop'])
#    promptHist.Add(histos['ww_incl'])
#    promptHist.Add(histos['wz'])
#    promptHist.Add(histos['zz_incl'])
#    promptHist.Add(histos['ztautau'])

    dibosonsHist = histos['ww_incl'].Clone('dibosons')
    dibosonsHist.Add(histos['wz'])
    dibosonsHist.Add(histos['zz_incl'])

    zdyHist = histos['dy50'].Clone('zdy')
    zdyHist.Add(histos['dy100to200'])
    zdyHist.Add(histos['dy200to400'])
    zdyHist.Add(histos['dy400to500'])
    zdyHist.Add(histos['dy500to700'])
    zdyHist.Add(histos['dy700to800'])
    zdyHist.Add(histos['dy800to1000'])
    zdyHist.Add(histos['dy1000to1500'])


#    zprime = histos['zssm1000'].Clone('zprime')

    # simulate stacking
#    promptHist.Add(jetsHist)
#    zdyHist.Add(jetsHist)
#    zdyHist.Add(promptHist)

    fexport = ROOT.TFile('histos_export%s.root' % ('_cumulative' if cumulative else ''), 'RECREATE')
    dataHist.Write()
    zdyHist.Write()
    promptHist.Write()
    dibosonsHist.Write()
    jetsHist.Write()
#    zprime.Write()
    fexport.Close()

    if not cumulative:
        print 'our Z peak events in data:', dataHist.Integral(dataHist.FindBin(60),dataHist.FindBin(120)-1)
        print 'our Z peak events in MC:  ', zdyHist.Integral(zdyHist.FindBin(60),zdyHist.FindBin(120)-1)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

## Now e-mu.
#histos = {}
#for sample in samples:
#    f = ROOT.TFile(emu_mc_fn_base % sample.name)
#    h = f.Histosstartup.Get(emu_histogram).Clone()
#    h.Rebin(emu_rebin_factor)
#    h.Scale(sample.partial_weight * emu_scale)
#    histos[sample.name] = h
#
#f = ROOT.TFile(emu_data_fn)
#dataHist = f.Histosstartup.Get(emu_histogram).Clone('dataHist')
#dataHist.Rebin(emu_rebin_factor)
#
#jetsHist = histos['inclmu15'].Clone('jetsHist')
#jetsHist.Add(histos['wjets'])
#
#promptHist = histos['ttbar'].Clone('promptHist')
#for x in ['tW', 'tbarW', 'ww', 'wz', 'zz', 'ztautau', 'zmumu', 'dy120', 'dy200', 'dy500', 'dy800']:
#    promptHist.Add(histos[x])
#
## Simulate stack.
#promptHist.Add(jetsHist)
#                
#fexport = ROOT.TFile('histos_export_emu.root', 'RECREATE')
#dataHist.Write()
#promptHist.Write()
#jetsHist.Write()
#fexport.Close()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

# In the past the HEEP histograms came pre-added, now they don't.
if not add_heep:
    sys.exit(0)

for cumulative in (False, True):
    f = ROOT.TFile(heep_fn)
    
    jetsHist = f.Get('qcdHistEBEB').Clone('jetsHist')
    jetsHist.Add(f.Get('qcdHistEBEE'))
    jetsHist.Add(f.Get('wjetHistEBEB'))
    jetsHist.Add(f.Get('wjetHistEBEE'))
    jetsHist.Add(f.Get('phoJetHistEBEB'))
    jetsHist.Add(f.Get('phoJetHistEBEE'))

    promptHist = f.Get('ttbarHistEBEB').Clone('promptHist')
    promptHist.Add(f.Get('ttbarHistEBEE'))

    zdyHist = f.Get('zeeHistEBEB').Clone('zdyHist')
    zdyHist.Add(f.Get('zeeHistEBEE'))

    dataHist = f.Get('dataHistEBEB').Clone('dataHist')
    dataHist.Add(f.Get('dataHistEBEE'))

    if heep_rebin_roundto is not None and not cumulative:
        for h in (jetsHist, promptHist, zdyHist, dataHist):
            variable_rebin_and_normalize(h, None, heep_rebin_roundto)
    elif heep_rebin_factor > 1:
        for h in (jetsHist, promptHist, zdyHist, dataHist):
            h.Rebin(heep_rebin_factor)
    
    if cumulative:
        jetsHist   = cumulative_histogram(jetsHist)
        promptHist = cumulative_histogram(promptHist)
        zdyHist    = cumulative_histogram(zdyHist)
        dataHist   = cumulative_histogram(dataHist)
        for h in [jetsHist, promptHist, zdyHist, dataHist]:
            h.SetName(h.GetName().replace('_cumulative_ge', ''))

    # simulate stacking
    promptHist.Add(jetsHist)
    zdyHist.Add(jetsHist)
    zdyHist.Add(promptHist)

    fexport = ROOT.TFile('histos_export_heep%s.root' % ('_cumulative' if cumulative else ''), 'RECREATE')
    dataHist.Write()
    zdyHist.Write()
    promptHist.Write()
    jetsHist.Write()
    fexport.Close()

    if not cumulative:
        print 'HEEP Z peak events in data:', dataHist.Integral(dataHist.FindBin(60), dataHist.FindBin(120)-1)
        print 'HEEP Z peak events in MC:  ', zdyHist .Integral(zdyHist .FindBin(60), zdyHist .FindBin(120)-1)
