#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import set_zp2mu_style, plot_saver, ROOT
set_zp2mu_style()

mumu_fn, ee_fn, = sys.argv[1:3]
mumu = ROOT.TFile(mumu_fn)
ee = ROOT.TFile(ee_fn)

mumu_data  = mumu.Get('VBTFMuonsPlusMuonsMinus_data')
mumu_sumMC = mumu.Get('VBTFMuonsPlusMuonsMinus_sumMC')

ee_data = ee.Get('dataHist')
ee_MC = ee.Get('zeeHist'), ee.Get('ttbarHist'), ee.Get('qcdHist')

binning = 5
eenb = ee_data.GetNbinsX()
assert((ee_data.GetBinLowEdge(eenb + 1) - ee_data.GetBinLowEdge(1))/eenb == binning)
for eemc in ee_MC:
    assert(eenb == eemc.GetNbinsX())
    assert((eemc.GetBinLowEdge(eenb + 1) - eemc.GetBinLowEdge(1))/eenb == binning)
mumunb = mumu_data.GetNbinsX()
assert(mumunb == mumu_sumMC.GetNbinsX())
assert((mumu_data.GetBinLowEdge(mumunb + 1) - mumu_data.GetBinLowEdge(1))/mumunb == binning)
assert((mumu_sumMC.GetBinLowEdge(mumunb + 1) - mumu_sumMC.GetBinLowEdge(1))/mumunb == binning)

total_data  = ROOT.TH1F('total_data', '', mumunb, mumu_data.GetXaxis().GetXmin(), mumu_data.GetXaxis().GetXmax())
total_sumMC = ROOT.TH1F('total_sumMC', '', mumunb, mumu_data.GetXaxis().GetXmin(), mumu_data.GetXaxis().GetXmax())

for mumu_bin in xrange(1, mumunb+2):
    mumu_bin_value = mumu_data.GetBinLowEdge(mumu_bin) + binning/10.
    ee_bin = ee_data.FindBin(mumu_bin_value)
    total_data.SetBinContent(mumu_bin, mumu_data.GetBinContent(mumu_bin) + ee_data.GetBinContent(ee_bin))
    total_data.SetBinError(mumu_bin, (mumu_data.GetBinError(mumu_bin)**2 + ee_data.GetBinError(ee_bin)**2)**0.5)

    mc = mumu_sumMC.GetBinContent(mumu_bin)
    mcvar = mumu_sumMC.GetBinContent(mumu_bin)**2

    for eemc in ee_MC:
        ee_bin = eemc.FindBin(mumu_bin_value)
        mc += eemc.GetBinContent(ee_bin)
        mcvar += eemc.GetBinError(ee_bin)**2

    total_sumMC.SetBinContent(mumu_bin, mc)
    total_sumMC.SetBinError(mumu_bin, mcvar**0.5)

for h in (total_data, total_sumMC):
    h.SetTitle('')
    h.GetXaxis().SetRangeUser(40, 400)
    h.GetYaxis().SetRangeUser(1e-2, 4e3)
    h.GetXaxis().SetTitle('dilepton (ee + #mu#mu) invariant mass (GeV)')
    h.GetYaxis().SetTitle('Events/5 GeV')
    h.SetStats(0)
    
ps = plot_saver('plots/datamc/combination')
total_data.Draw()
ps.save('data')
total_sumMC.Draw()
ps.save('sumMC')

total_data.SetMarkerStyle(20)
total_data.SetMarkerSize(0.5)
total_sumMC.SetLineColor(ROOT.kBlue)
total_data.Draw('e1')
total_sumMC.Draw('same hist')
l = ROOT.TLegend(0.62, 0.71, 0.87, 0.87)
l.SetFillColor(0)
l.AddEntry(total_data, 'Data', 'LP')
l.AddEntry(total_sumMC, 'SM background', 'L')
l.Draw('same')
t1 = ROOT.TLatex(0.4, 0.93, '#sqrt{s} = 7 TeV,  #int L dt = 8.5 pb^{-1}')
t2 = ROOT.TLatex(0.1, 0.93, 'CMS preliminary')

for t in t1, t2:
    t.SetTextSize(0.0375)
    t.SetNDC()

t1.Draw()
#t2.Draw()

ps.save('comparison')

for h in (total_data, total_sumMC):
    h.GetYaxis().SetRangeUser(0, 2800)
    h.GetYaxis().SetLabelSize(0.03)
total_data.Draw('e1')
total_sumMC.Draw('same hist')
l.Draw('same')
t1.Draw()
#t2.Draw()
ps.save('comparison_lin')

def cumulative(h):
    nb = h.GetNbinsX()
    hc = ROOT.TH1F(h.GetName() + '_cumulative', '', nb, h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    for i in xrange(nb+1, 0, -1):
        prev = 0 if i == nb+1 else hc.GetBinContent(i+1)
        hc.SetBinContent(i, h.GetBinContent(i) + prev)
    return hc

total_data_cum = cumulative(total_data)
total_sumMC_cum = cumulative(total_sumMC)
total_sumMC_cum.SetLineColor(ROOT.kBlue)
for h in (total_data_cum, total_sumMC_cum):
    h.SetTitle('')
    h.GetXaxis().SetRangeUser(40, 400)
    h.GetYaxis().SetRangeUser(0.3, 1e4)
    h.GetXaxis().SetTitle('dilepton (ee + #mu#mu) invariant mass (GeV)')
    h.GetYaxis().SetTitle('Events < X')
    h.SetStats(0)
total_data_cum.Draw()
total_sumMC_cum.Draw('same')
l = ROOT.TLegend(0.62, 0.71, 0.87, 0.87)
l.SetFillColor(0)
l.AddEntry(total_data_cum, 'Data', 'L')
l.AddEntry(total_sumMC_cum, 'SM background', 'L')
l.Draw('same')
t1.Draw()
#t2.Draw()
ps.save('cumulative')
for h in (total_data_cum, total_sumMC_cum):
    h.GetYaxis().SetRangeUser(0, 6000)
    h.GetYaxis().SetLabelSize(0.03)
total_data_cum.Draw()
total_sumMC_cum.Draw('same hist')
l.Draw('same')
t1.Draw()
#t2.Draw()
ps.save('cumulative_lin')

            
