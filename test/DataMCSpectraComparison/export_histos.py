#!/usr/bin/env python

# Script to export our histograms into the same format that HEEP uses
# so they can be drawn with a common script (plot_for_paper.C).

import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *

ROOT.TH1.AddDirectory(False)

mumu_mc_fn_base = 'mc/ana_datamc_%s.root'
mumu_data_fn = 'data/Run2011MuonsOnly/ana_datamc_data.root'
mumu_rebin_factor = 5
mumu_histogram = 'DimuonMassVertexConstrained'
mumu_scale = 4914*1.0969257655
emu_mc_fn_base = 'mc/ana_datamc_%s.root'
emu_data_fn = 'data/Run2011/ana_datamc_data.root'
emu_rebin_factor = 20
emu_histogram = 'DileptonMass'
emu_scale = 4641
add_heep = True
heep_fn = 'heep_massHists2011Jan9.root'
heep_rebin_factor = None

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

# First mu+mu-, differential and cumulative.
for cumulative in (False, True):
    histos = {}
    for sample in samples:
        f = ROOT.TFile(mumu_mc_fn_base % sample.name)
        h = f.OurNewMuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone()
        h.Rebin(mumu_rebin_factor)
        h.Scale(sample.partial_weight * mumu_scale)
        if cumulative:
            h = cumulative_histogram(h)
            h.SetName(h.GetName().replace('_cumulative_ge', ''))
        histos[sample.name] = h

    f = ROOT.TFile(mumu_data_fn)
    dataHist = f.OurNewMuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone('dataHist')
    dataHist.Rebin(mumu_rebin_factor)
    if cumulative:
        dataHist = cumulative_histogram(dataHist)
        dataHist.SetName('dataHist')
    
    jetsHist = histos['inclmu15'].Clone('jetsHist')
    jetsHist.Add(histos['wjets'])

    promptHist = histos['ttbar'].Clone('promptHist')
    promptHist.Add(histos['tW'])
    promptHist.Add(histos['tbarW'])
    promptHist.Add(histos['ww'])
    promptHist.Add(histos['wz'])
    promptHist.Add(histos['zz'])
    promptHist.Add(histos['ztautau'])
    
    zdyHist = histos['zmumu'].Clone('zdyHist')
    zdyHist.Add(histos['dy120'])
    zdyHist.Add(histos['dy200'])
    zdyHist.Add(histos['dy500'])
    zdyHist.Add(histos['dy800'])
    zdyHist.Add(histos['dy1000'])

    zprime = histos['zssm1000'].Clone('zprime')
    
    # simulate stacking
    promptHist.Add(jetsHist)
    zdyHist.Add(jetsHist)
    zdyHist.Add(promptHist)

    fexport = ROOT.TFile('histos_export%s.root' % ('_cumulative' if cumulative else ''), 'RECREATE')
    dataHist.Write()
    zdyHist.Write()
    promptHist.Write()
    jetsHist.Write()
    zprime.Write()
    fexport.Close()

    if not cumulative:
        print 'our Z peak events in data:', dataHist.Integral(dataHist.FindBin(60),dataHist.FindBin(120)-1)
        print 'our Z peak events in MC:  ', zdyHist.Integral(zdyHist.FindBin(60),zdyHist.FindBin(120)-1)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

# Now e-mu.
histos = {}
for sample in samples:
    f = ROOT.TFile(emu_mc_fn_base % sample.name)
    h = f.OurNewMuonsElectronsOppSignHistos.Get(emu_histogram).Clone()
    h.Rebin(emu_rebin_factor)
    h.Scale(sample.partial_weight * emu_scale)
    histos[sample.name] = h

f = ROOT.TFile(emu_data_fn)
dataHist = f.OurNewMuonsElectronsOppSignHistos.Get(emu_histogram).Clone('dataHist')
dataHist.Rebin(emu_rebin_factor)

jetsHist = histos['inclmu15'].Clone('jetsHist')
jetsHist.Add(histos['wjets'])

promptHist = histos['ttbar'].Clone('promptHist')
for x in ['tW', 'tbarW', 'ww', 'wz', 'zz', 'ztautau', 'zmumu', 'dy120', 'dy200', 'dy500', 'dy800']:
    promptHist.Add(histos[x])

# Simulate stack.
promptHist.Add(jetsHist)
                
fexport = ROOT.TFile('histos_export_emu.root', 'RECREATE')
dataHist.Write()
promptHist.Write()
jetsHist.Write()
fexport.Close()

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

    if heep_rebin_factor > 1:
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
