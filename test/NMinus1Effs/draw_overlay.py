#!/usr/bin/env python

# Draws nminus1 histos from two files overlayed.
# Both data and MC, M=60-120 GeV bin.

from pprint import pprint
import sys, os
import ROOT
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
#from roottools import *
from ROOT import *
set_zp2mu_style()

psn = 'plots/nminus1effs_overlay'
ps = plot_saver(psn, size=(900,600), log=False, pdf=True)
ps.c.SetBottomMargin(0.2)

infile1 = ROOT.TFile("whargl_oldid.root","open")
infile2 = ROOT.TFile("whargl_newid.root","open")

f1_data = infile1.arp0
f2_data = infile2.arp0
f1_mc   = infile1.arp1
f2_mc   = infile2.arp1

f2_data.SetTitle('')
f2_data.SetMinimum(0.7)
f2_data.SetMaximum(1.02)
f2_data.SetMarkerStyle(24)
f2_data.SetMarkerColor(kBlue)
f2_data.SetLineColor(kBlue)
f2_data.SetLineStyle(2)
f2_data.Draw("ape")

f1_data.SetMarkerStyle(20)
f1_data.SetMarkerColor(kBlack)
f1_data.SetLineColor(kBlack)
f1_data.SetLineStyle(1)
f1_data.Draw("pe, same")

f1_mc.SetLineColor(kBlack)
f1_mc.SetLineStyle(1)
f1_mc.Draw("e,same")
f2_mc.SetLineColor(kBlue)
f2_mc.SetLineStyle(2)
f2_mc.Draw("e,same")

leg = TLegend()
leg.SetFillColor(0)
leg.SetTextAlign(12)
leg.SetTextSize(0.04)
leg.SetX1NDC(0.72)
leg.SetY1NDC(0.750)
leg.SetX2NDC(0.82)
leg.SetY2NDC(0.85)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

leg.SetX1NDC(0.55)
leg.SetY1NDC(0.32)
leg.SetX2NDC(0.75)
leg.SetY2NDC(0.52)

leg.AddEntry(f1_data, "data, old muon id", "lp")
leg.AddEntry(f2_data, "data, new muon id", "lp")
leg.AddEntry(f1_mc,   "MC, old muon id", "lp")
leg.AddEntry(f2_mc,   "MC, new muon id","l")
leg.Draw("same")

tex = TLatex()
tex.SetTextColor(kBlack)
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetTextAlign(22)
tex.SetNDC()
tex.DrawLatex(0.70,0.57,"60 < M < 120 GeV")

ps.save("draw_overlay")
