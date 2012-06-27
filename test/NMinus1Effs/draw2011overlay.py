#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

from pprint import pprint
import sys, os
import ROOT
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
#from roottools import *
from ROOT import *
set_zp2mu_style()

infile = ROOT.TFile("whargl.root","open")
infile1 = ROOT.TFile("data2011out.root","open")
can = TCanvas("can","can",600,600)


gre1 = infile.arp0
gre2 = infile.arp1

gre1.SetMinimum(0.7)
#gre1.SetMaximum(1.1)
can.SetBottomMargin(0.2)
gre1.Draw("ape")
gre2.Draw("pe,same")
can.Update()
gre1.Draw("ape")
gre2.Draw("pe,same")


gre3 = gre1.Clone("gre3")
gre3a=infile1.gre
print gre3a.GetN(),gre3.GetN()
for ip in range(1,gre1.GetN()+1):
	ip1 = ip-1

	nv = gre3a.GetY()[ip1] 
	ne = gre3a.GetErrorY(ip1) 

	xv = gre3.GetX()[ip]	
	xel = gre3.GetErrorXlow(ip)	
	xeh = gre3.GetErrorXhigh(ip)	
	
	gre3.SetPoint(ip,xv,nv)
	gre3.SetPointError(ip,xel,xeh,0.5*ne,0.5*ne)

gre3.SetLineColor(kBlue)
gre3.SetMarkerColor(kBlue)
gre3.SetMarkerStyle(23)

gre1.SetTitle('')
gre1.Draw("ape")
gre2.Draw("pe,same")
gre3.Draw("pe,same")

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

leg.SetX1NDC(0.20)
leg.SetY1NDC(0.22)
leg.SetX2NDC(0.28)
leg.SetY2NDC(0.35)

leg.AddEntry(gre1,"data 2012","lp")
leg.AddEntry(gre3,"data 2011","lp")
leg.AddEntry(gre2,"MC 2012","lp")
#leg.AddEntry(gre2,"MC 2012","l")
leg.Draw("same")

tex = TLatex()
tex.SetTextColor(kBlack)
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetTextAlign(22)
tex.SetNDC()

tex.DrawLatex(0.32,0.38,"60 < M < 120")

can.Update()

