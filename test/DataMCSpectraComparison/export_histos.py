#!/usr/bin/env python

# Script to export our histograms into the same format that HEEP uses
# so they can be drawn with a common script (plot_for_paper.C).

import sys
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from samples import *

ROOT.TH1.AddDirectory(False)

data_fn_base = '/uscms/home/tucker/nobackup/ana_datamc_V00-10-26'
mc_fn_base = '/uscms/home/tucker/nobackup/zp2mu_ana_datamc_mc/V00-10-26'
mumu_mc_fn_base = os.path.join(mc_fn_base, 'ana_datamc_%s.root')
mumu_data_fn = os.path.join(data_fn_base, 'ana_datamc_Run2011AMuonsOnly/ana_datamc_data.root')
mumu_rebin_factor = 5
mumu_histogram = 'DimuonMassVertexConstrained'
mumu_scale = 1125.74317154*1.1603697643990654
emu_mc_fn_base = os.path.join(mc_fn_base, 'ana_datamc_%s.root')
emu_data_fn = os.path.join(data_fn_base, 'ana_datamc_Run2011A/ana_datamc_data.root')
emu_rebin_factor = 20
emu_histogram = 'DileptonMass'
emu_scale = 1086.9
add_heep = True
heep_fn = 'massHist1100pb.root'

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
        histos[sample.name] = h

    f = ROOT.TFile(mumu_data_fn)
    dataHist = f.OurNewMuonsPlusMuonsMinusHistos.Get(mumu_histogram).Clone('dataHist')
    dataHist.Rebin(mumu_rebin_factor)
    if cumulative:
        dataHist = cumulative_histogram(dataHist)
        dataHist.SetName('dataHist')
    
    qcdHist = histos['inclmu15'].Clone('qcdHist')
    qcdHist.Add(histos['wjets'])

    ttbarHist = histos['ttbar'].Clone('ttbarHist')
    ttbarHist.Add(histos['singletop_tW'])
    ttbarHist.Add(histos['ww'])
    ttbarHist.Add(histos['wz'])
    ttbarHist.Add(histos['zz'])
    ttbarHist.Add(histos['ztautau'])
    
    zeeHist = histos['zmumu'].Clone('zeeHist')
    zeeHist.Add(histos['dy200'])
    zeeHist.Add(histos['dy500'])
    zeeHist.Add(histos['dy800'])
    zeeHist.Add(histos['dy1000'])

    zprime = histos['zssm1000'].Clone('zprime')
    
    # simulate stacking
    ttbarHist.Add(qcdHist)
    zeeHist.Add(qcdHist)
    zeeHist.Add(ttbarHist)

    fexport = ROOT.TFile('histos_export%s.root' % ('_cumulative' if cumulative else ''), 'RECREATE')
    dataHist.Write()
    zeeHist.Write()
    ttbarHist.Write()
    qcdHist.Write()
    zprime.Write()
    fexport.Close()

    if not cumulative:
        print 'our Z peak events in data:', dataHist.Integral(dataHist.FindBin(60),dataHist.FindBin(120)-1)
        print 'our Z peak events in MC:  ', zeeHist.Integral(zeeHist.FindBin(60),zeeHist.FindBin(120)-1)

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
for x in ['zmumu', 'singletop_tW', 'ww', 'wz', 'zz', 'ztautau', 'dy200', 'dy500', 'dy800']:
    promptHist.Add(histos[x])
promptHist.Add(jetsHist)

zeeHist = histos['zmumu'].Clone('zeeHist')
zeeHist.Add(histos['dy200'])
zeeHist.Add(histos['dy500'])
zeeHist.Add(histos['dy800'])
zeeHist.Add(histos['dy1000'])
                
fexport = ROOT.TFile('histos_export_emu.root', 'RECREATE')
dataHist.Write()
jetsHist.Write()
promptHist.Write()
fexport.Close()

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

    ttbarHist = f.Get('ttbarHistEBEB').Clone('ttbarHist')
    ttbarHist.Add(f.Get('ttbarHistEBEE'))

    zeeHist = f.Get('zeeHistEBEB').Clone('zeeHist')
    zeeHist.Add(f.Get('zeeHistEBEE'))

    dataHist = f.Get('dataHistEBEB').Clone('dataHist')
    dataHist.Add(f.Get('dataHistEBEE'))
    
    if cumulative:
        jetsHist  = cumulative_histogram(jetsHist)
        ttbarHist = cumulative_histogram(ttbarHist)
        zeeHist   = cumulative_histogram(zeeHist)
        dataHist  = cumulative_histogram(dataHist)

    # simulate stacking
    ttbarHist.Add(jetsHist)
    zeeHist.Add(jetsHist)
    zeeHist.Add(ttbarHist)

    fexport = ROOT.TFile('histos_export_heep%s.root' % ('_cumulative' if cumulative else ''), 'RECREATE')
    dataHist.Write()
    zeeHist.Write()
    ttbarHist.Write()
    jetsHist.Write()
    fexport.Close()

    if not cumulative:
        print 'HEEP Z peak events in data:', dataHist.Integral(dataHist.FindBin(60),dataHist.FindBin(120)-1)
        print 'HEEP Z peak events in MC:  ', zeeHist.Integral(zeeHist.FindBin(60),zeeHist.FindBin(120)-1)
