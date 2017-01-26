#!/usr/bin/env python

# Draws nminus1 histos from two files overlayed.
# Both data and MC, M=60-120 GeV bin.

from pprint import pprint
import sys, os
import ROOT
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
#from roottools import *
from ROOT import *
#set_zp2mu_style()
gStyle.SetOptStat("")

titleoffset = 1.1
textfont = 40
textsize = 0.05
#newStyle = ROOT.TStyle("EXAMPLE", "Example style");
gStyle.SetPadTopMargin(0.005)
gStyle.SetPadRightMargin(0.005)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetTitleXOffset(titleoffset)
gStyle.SetTitleYOffset(titleoffset)
gStyle.SetTextFont(textfont)
gStyle.SetTextSize(textsize)
#gROOT.SetStyle("EXAMPLE")
#gROOT.ForceStyle()

infile1 = ROOT.TFile("dimuon_histos_differential_madgraph.root","open")
mad = infile1.Get('zdy').Clone()
infile2 = ROOT.TFile("dimuon_histos_differential_powheg.root","open")
pow = infile2.Get('zdy').Clone()
#from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import dy50, dy50to120
#lumi = 48
#print dy50.name, dy50.partial_weight*lumi
#print dy50to120.name, dy50to120.partial_weight*lumi
#mad.Scale(dy50.partial_weight * lumi)
#pow.Scale(dy50to120.partial_weight * lumi)
mad.SetTitle(" ")
mad.SetLineColor(8)
mad.SetLineWidth(2)
pow.SetLineColor(1)
pow.SetLineWidth(2)

fair_share = 0.7


c = ROOT.TCanvas( 'c', '',800, 800)
c.SetFrameFillStyle(0)
c.cd()
pad1 = ROOT.TPad("pad1","pad1",0,1-fair_share,1,1)
pad1.Draw()
pad1.cd()
pad1.SetLogy()
pad1.SetBottomMargin(0)
mad.GetXaxis().SetRangeUser(40,1300)
pow.GetXaxis().SetRangeUser(40,1300)
#mad.GetXaxis().SetRangeUser(260,350)
#pow.GetXaxis().SetRangeUser(260,350)
mad.Draw("e")
pow.Draw("e same ")

leg = ROOT.TLegend( 0.55, 0.75, 0.85, 0.9 )
leg.SetFillColor(kWhite)
#leg.SetLineColor(kWhite)
leg.SetTextSize(0.03)
leg.AddEntry(mad, " DY madgraph", "epl")
leg.AddEntry(pow, " DY powheg", "epl")
leg.Draw('same')
c.cd()
pad2 = ROOT.TPad("pad2","pad2",0,0,1,1-fair_share)
pad2.SetTopMargin(0)
pad2.SetGrid()
pad2.Draw()
pad2.cd()
ratio = mad.Clone()
ratio.Sumw2()
ratio.Divide(pow.Clone())
ratio.Draw("pe")
ratio.GetYaxis().SetTitle("madgraph/powheg")
ratio.GetYaxis().SetNdivisions(505)
#// TAttAxis::SetLabelSize says:  The size is expressed in per cent of the pad width
#// same do TAttAxis::SetTitleSize, TAttAxis::SetLabelOffset and TAttAxis::SetTitleOffset...
ratio.GetXaxis().SetLabelSize(textsize * fair_share / (1 - fair_share))
ratio.GetXaxis().SetTitleSize(textsize * fair_share / (1 - fair_share))
ratio.GetXaxis().SetLabelSize(textsize * fair_share / (1 - fair_share))
ratio.GetYaxis().SetLabelSize(textsize * fair_share / (1 - fair_share))
ratio.GetYaxis().SetTitleSize(textsize * fair_share / (1 - fair_share))
ratio.GetYaxis().SetTitleOffset(titleoffset * fair_share / (1 - fair_share))
ratio.GetXaxis().SetLabelSize(0.08)
ratio.GetYaxis().SetLabelSize(0.08)
ratio.GetYaxis().SetTitleSize(0.08)
ratio.GetYaxis().SetTitleOffset(0.5)
ratio.GetYaxis().SetRangeUser(0.,2.)
ratio.SetLineColor(4)

c.Print('pow_mad.pdf')


