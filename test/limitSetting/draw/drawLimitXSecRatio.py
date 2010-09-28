#!/usr/bin/env py

import ROOT, commands, math
from ROOT import *

from ufPlotTools import *
setTDRStyle()


sFileOfInterest	= "out400b"
sDataFile		=	"dimuons_merged_2900InvNb_nocosmics_new.root"

fileData = TFile("data/%s"%sDataFile,"open")

#can = TCanvas("can","can", 1000,600)
can = TCanvas("can","can", 400,400)
listToSave = []
can.Draw()
'''
treeData = fileData.Get("tree")
zhist = TH1F("zhist", "", 20, 60, 120)
zhist.Sumw2()
zhist.SetXTitle("#mu#mu invariant mass")
zhist.SetYTitle("entries / %2.1f GeV/c^{2}"%zhist.GetXaxis().GetBinWidth(1))
treeData.Project("zhist","recoCandMass","reco1.charge != reco2.charge & recoCandMass>60 & recoCandMass<=120")

numZinSample = zhist.Integral(1,20)
#print numZinSample
#print zhist.GetEntries()
zhist.SetTitle("estimated number of Zs: %4.f"%numZinSample)
zhist.Draw()
'''
def GetZHistogram(fileData):
	treeData = fileData.Get("tree")
	
	zhist =TH1F("zhist", "", 20, 60, 120)
	zhist.Sumw2()
	zhist.SetXTitle("#mu#mu invariant mass")
	zhist.SetYTitle("entries / %2.1f GeV/c^{2}"%zhist.GetXaxis().GetBinWidth(1))
	treeData.Project("zhist","recoCandMass","reco1.charge != reco2.charge & recoCandMass>60 & recoCandMass<=120")
	
	numZinSample = zhist.Integral(1,20)
	#print numZinSample
	#print zhist.GetEntries()
	zhist.SetTitle("estimated number of Zs: %4.f"%numZinSample)
	return zhist
	
#can.SetWindowSize(400,400)
zhist = GetZHistogram(fileData)
zhist.SetName("zpole")
#zhist.SetXTitle("di-muon mass [GeV/c^{2}]")
zhist.Draw("")
can.Update()
can.SaveAs("gif/zpole.gif")
#can.SetWindowSize(1000,600)

#listToSave.append(zhist)

print zhist.GetEntries()

import EffRatio

greEff = TGraphErrors()

zprimemasses = [500,750,1000,1250,1500,1750]

for mass in zprimemasses:
	ipoint = greEff.GetN()
	greEff.SetPoint(ipoint		,mass	,EffRatio.ZPrimeEff(mass))
	greEff.SetPointError(ipoint	,0		, EffRatio.ZPrimeEffErrBinomial(mass))
	

greEff.SetMinimum(0)
greEff.SetMaximum(1)
greEff.GetXaxis().SetTitle("z-prime mass")
greEff.GetYaxis().SetTitle("z-prime efficiency")

greEff.Draw("ape")
greEff.SetName("eff")
listToSave.append(greEff)

greEffRatio = TGraphErrors()
for mass in zprimemasses:
	ipoint = greEffRatio.GetN()
	greEffRatio.SetPoint(ipoint		,mass	,EffRatio.GetEffRatio(mass))

greEffRatio.SetMinimum(0)
greEffRatio.SetMaximum(1)
greEffRatio.GetXaxis().SetTitle("z-prime mass")
#greEffRatio.GetYaxis().SetTitle("#epsilon_{z}/#epsilon_{z^{{}^{#acute}}}")#")#_{z-prime}")
greEffRatio.GetYaxis().SetTitle("#epsilon_{z}/#epsilon_{z-prime}")
greEffRatio.Draw("ape")

greEffRatio.SetName("effRatio")
listToSave.append(greEffRatio)

greEffRatioExtrap = TGraphErrors()
greEffRatioExtrap.SetMinimum(0)
greEffRatioExtrap.SetMaximum(2)
greEffRatioExtrap.GetXaxis().SetTitle("z-prime mass")
greEffRatioExtrap.GetYaxis().SetTitle("extrapolated #epsilon_{z}/#epsilon_{z-prime}")
for mass in range(100,1750):
	ipoint = greEffRatioExtrap.GetN()
	greEffRatioExtrap.SetPoint(ipoint		,mass	,EffRatio.GetEffRatioExtrap(mass))
	

greEffRatioExtrap.SetName("greEffRatioExtrap")
listToSave.append(greEffRatioExtrap)

#
#
#
#
def GetXSecRatioQuick(lambdaS, numZ, effRatio ):
#	return lambdaS*effZ/(effZPrime*Nz)
	return lambdaS*effRatio/(1.*numZ)

infileLimit	= TFile("out/limitPlot_%s.root"%sFileOfInterest, "open")

limitPlot = infileLimit.Get("limitVsMass")
#limitPlot.Draw("ap")
#can.Update()

limitRatio = limitPlot.Clone("limitRatio")
limitRatio.SetMarkerColor(kRed)
#limitRatio.Draw("ape")
#can.Update()

def GetXSecRatioFromTGraph(graph,ipoint,zhist):
	xval = limitRatio.GetX()[ipoint]
	yval = limitRatio.GetY()[ipoint]
	lambdaS = yval
	mass = xval
	effRatio = EffRatio.GetEffRatioExtrap(mass)
	numZ = zhist.GetEntries() 
#	numZ = 30000.
#	numZ = 10000
	xsecRat = GetXSecRatioQuick(lambdaS, numZ, effRatio )
	print "%d\t%f\t%f"%(ipoint,xval,xsecRat)
	graph.SetPoint(ipoint,xval,xsecRat)
#	return GetXSecRatioQuick(lambdaS, numZ, effRatio )

numPoints = limitRatio.GetN()
for ipoint in range(0,numPoints):
	GetXSecRatioFromTGraph(limitRatio,ipoint,zhist)

limitRatio.GetXaxis().SetTitle("z-prime mass")
limitRatio.GetYaxis().SetTitle("limit on #sigma(z-prime)/#sigma(z)")
limitRatio.SetTitle("limit on the Z^{_{/}}/Z cross-section ratio")

limitRatio.SetMaximum(0.0016)
limitRatio.SetMarkerSize(2)
limitRatio.Draw("ap")
can.Update()

can.SaveAs("gif/limitRatio.gif")

#	yval = 0
#	limitRatio.GetPoint(ipoint,xval,yval)


#	greEff.SetPointError(ipoint	,0		, EffRatio.ZPrimeEffErrBinomial(mass))

# 
# expected Z' plot
# 
# 

'''

intLumi = 30.0

greExpecteLambdaS_s6 = TGraphErrors()
greExpecteLambdaS_s6.SetPoint(0,1000,intLumi*55e-3)
greExpecteLambdaS_s6.SetPoint(1,1500,intLumi*4.6e-3)
greExpecteLambdaS_s6.SetPoint(2,2000,intLumi*0.46e-3)





greExpecteLambdaS_s10 = TGraphErrors()
greExpecteLambdaS_s10.SetPoint(0,1000,intLumi*220e-3)
greExpecteLambdaS_s10.SetPoint(1,1500,intLumi*22e-3)
greExpecteLambdaS_s10.SetPoint(2,2000,intLumi*6.4e-3)



greExpecteLambdaS_s6.SetMinimum(0)
greExpecteLambdaS_s6.SetMaximum(3)
greExpecteLambdaS_s10.SetMarkerColor(kRed)

greExpecteLambdaS_s6.Draw("ape")
greExpecteLambdaS_s10.Draw("p,same")

numExpectedZ = 1686 * intLumi 

print numExpectedZ

greExpecteLambdaS_s10_lim = TGraphErrors()
mass = 1000
greExpecteLambdaS_s10_lim.SetPoint(0,mass,intLumi*220e-3/numExpectedZ*EffRatio.GetEffRatioExtrap(mass))
mass = 1500
greExpecteLambdaS_s10_lim.SetPoint(1,mass,intLumi*22e-3/numExpectedZ*EffRatio.GetEffRatioExtrap(mass))
mass = 2000
greExpecteLambdaS_s10_lim.SetPoint(2,mass,intLumi*6.4e-3/numExpectedZ*EffRatio.GetEffRatioExtrap(mass))

greExpecteLambdaS_s10_lim.Draw("ape")

limitRatio.Draw("ape")
'''
fit = TF1("fit", "expo", 100, 1800)
fit.SetParameters(-6.68568e0,-3.47780e-03)
fit.SetLineColor(kBlue)
fit.SetLineStyle(2)
fit.SetLineWidth(5)
fit.Draw("same")

can.Update()
can.SaveAs("gif/limitRatioWithTheory.gif")
#220, 33, 6.4

for plot in listToSave:
	plot.Draw("ape")
	can.Update()
	can.SaveAs("gif/%s.gif"%plot.GetName())


