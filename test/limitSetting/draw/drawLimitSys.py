#!/usr/bin/env py

import ROOT, commands, math
from ROOT import *

from ufPlotTools import *
setTDRStyle()

theDir = "out400b"
theDir = "sys400b"

files = commands.getoutput("ls -l %s | awk {'print $9'}"%theDir).rsplit("\n")
files.remove(files[0]) #remove the empty first part

def DumpFiles(files):
	for ifile in files:
		print ifile


def DumpFilesParse(files):
	for ifile in files:
		derp = ifile.rstrip(".root")
		print derp

#DumpFilesParse(files)

def GetInfoTuple(ifile):
	derp = ifile.rstrip(".root")
#	print derp
	derp = derp.split("_")
#	print derp	
	return derp
	
a = GetInfoTuple(files[0])
mass = a[1]
sysLambdab 		= int(a[3].strip("n"))
sysZPrimeWidth	= int(a[4].strip("n"))
#print sysLambdab

import ROOT
from ROOT import *


shits = []

def IsWanted(ifile,lambdaBVal,widthVal):
	info = GetInfoTuple(ifile)
	mass = info[1]
	sysLambdaB 		= int(info[3].strip("n"))
	sysZPrimeWidth	= int(info[4].strip("n"))

	if sysLambdaB != lambdaBVal: return False 
	if sysZPrimeWidth != widthVal:return False 
	return True


def GetFiles(files,lambdaBVal,widthVal):
	thefiles = []	
	for ifile in files:
		if not IsWanted(ifile,lambdaBVal,widthVal): continue
		thefiles.append(ifile)
	return thefiles

thefiles = GetFiles(files,10,-10)
#print thefiles

from ROOT import TFile,TGraphErrors,TF1
def GetGraphs(files):
	gres = []
	for ifile in files:
		infile = TFile("%s/%s"%(theDir,ifile),"open")
		newname = ifile.rstrip(".root")	
		print newname
		gre = infile.Get("pvalueVsLambdaS").Clone(newname)
		gre.SetTitle(newname);
		gres.append(gre)	
	return gres
		
#gres = GetGraphs(thefiles)
#print len(gres)

'''


for gre in gres:
	print gre.GetName()
	gre.Draw("ape")
	fit.SetParameters(-1,-1)
	gre.Fit("fit0","brq")
	can.Update()
	masspoint = int(gre.GetName().split("_")[1])
	limit = (math.log(0.05)-fit.GetParameter(0))/fit.GetParameter(1)

	print "%2d\t%4d\t%4.3f"%(ifit,masspoint,limit)
#	if masspoint <400: continue

	limitVsMass.SetPoint(ifit,masspoint,limit)
	ifit += 1	
'''

#limitVsMass.SetMinimum(0)
#limitVsMass.SetMaximum(7)
#limitVsMass.Draw("ape")
#can.Update()

ifit =0
limitVsMass = TGraphErrors()
fit = TF1("fit%d"%ifit,"expo",1,7)
fit.SetParameters(-1,-1)
can = TCanvas("can","can", 700,400)

def GetLimitPlot(gres,limitVsMass,fit):
	ifit = 0
	for gre in gres:
#		print gre.GetName()
		gre.Draw("ape")
		fit.SetParameters(-1,-1)
		gre.Fit("fit0","brq")
		can.Update()
		masspoint = int(gre.GetName().split("_")[1])
		limit = (math.log(0.05)-fit.GetParameter(0))/fit.GetParameter(1)

#		print "%2d\t%4d\t%4.3f"%(ifit,masspoint,limit)
#	if masspoint <400: continue
		ipoint = limitVsMass.GetN()
		limitVsMass.SetPoint(ipoint,masspoint,limit)
	

#GetLimitPlot(gres,limitVsMass,fit)

#limitVsMass.Draw("ape")	
#can.Update()


jitsWidth   = [-10,0,10]
jitsLambdaB = [-10,0,10]

limits = []

cols = [
	kBlack,	
	kRed,kBlue,kGreen,
	kYellow-2,kMagenta,kOrange,
	kAzure-2,kRed-2,kOrange+5
	]



icol = 0
for iwidth in jitsWidth:
	for ilambdaB in jitsLambdaB:
		name = "limit_%d_%d"%(ilambdaB,iwidth)
		limit = TGraphErrors()
		limit.SetName(name)	
		thefiles = GetFiles(files,ilambdaB,iwidth)
		gres = GetGraphs(thefiles)
#		print len(gres)
		
		limit.SetMarkerStyle(20)
		limit.SetMarkerSize(0.7)
		limit.SetMarkerColor(cols[icol])
		limit.SetFillColor(cols[icol])
		limit.SetLineColor(cols[icol])
		GetLimitPlot(gres,limit,fit)
		limits.append(limit)
		icol += 1
	
print "number of limits: %d"%len(limits)

	
	
limits[0].SetTitle("limit sytematics on #lambda_{B} and #Gamma")	
limits[0].GetXaxis().SetTitle("mass [GeV/c^{2}]")
limits[0].GetYaxis().SetTitle("#lambda_{s}")
limits[0].SetMinimum(0)
limits[0].SetMinimum(2)
limits[0].SetMaximum(8)
limits[0].Draw("ape")

icol =0
leg = TLegend()
leg.SetFillColor(0)
leg.SetTextAlign(12)
leg.SetX1NDC(0.5)
leg.SetY1NDC(0.4)
leg.SetX2NDC(0.8)
leg.SetY2NDC(0.85)

for limit in limits:
	limit.Draw("pe,same")
	sparse = limit.GetName().split("_")
#	print len(sparse)
#	print sparse[1]
#	print sparse[2]
	sentry = "#Delta#lambda_{S} = %+3d%-4s #Delta #Gamma = %+3d%s"%(int(sparse[1]),"%",int(sparse[2]),"%") 
#	print sentry
	leg.AddEntry(limit,sentry, "f")



leg.Draw("same")
can.Update()

can.SaveAs("gif/limitSys.gif")
can.SaveAs("gif/limitSys.pdf")


'''

print limits[0].GetName()

limits[1].SetMarkerColor(kRed)
limits[1].Draw("pe,same")
print limits[1].GetName()

limits[5].SetMarkerColor(kBlue)
limits[5].Draw("pe,same")
print limits[5].GetName()

#limits[14].SetMarkerColor(kGreen)
#limits[14].Draw("pe,same")
#print limits[14].GetName()
'''

'''

gres = []

for ifile in files:
#	path = "%s/%s"%(theDir,ifile)
#	print path
	infile = TFile("%s/%s"%(theDir,ifile),"open")
	newname = ifile.rstrip(".root")	
	print newname
	gre = infile.Get("pvalueVsLambdaS").Clone(newname)
	gre.SetTitle(newname);
	gres.append(gre)	
	

ifit = 0

fit = TF1("fit%d"%ifit,"expo",1,7)
fit.SetParameters(-1,-1)

can = TCanvas("can","can", 500,300)

gres.sort()

limitVsMass = TGraphErrors()

for gre in gres:
#	can.cd()
	gre.Draw("ape")
	fit.SetParameters(-1,-1)
	gre.Fit("fit0","brq")
#	gre.Fit("fit0","brq")
	can.Update()
#	can.SaveAs("gif/%s.gif"%gre.GetName())
#	print gre.GetName();
	masspoint = int(gre.GetName().split("_")[1])
	limit = (math.log(0.05)-fit.GetParameter(0))/fit.GetParameter(1)

	print "%2d\t%4d\t%4.3f"%(ifit,masspoint,limit)
#	print masspoint
	if masspoint <400: continue

	limitVsMass.SetPoint(ifit,masspoint,limit)
	ifit += 1	


limitVsMass.Draw("ape")

limitVsMass.GetXaxis().SetTitle("mass [GeV/c^{2}]")
limitVsMass.GetYaxis().SetTitle("#lambda_{s}")
limitVsMass.SetMinimum(0)

limitVsMass.Draw("ape")
can.Update()

line = TLine(limitVsMass.GetXaxis().GetXmin(),2.5,limitVsMass.GetXaxis().GetXmax(),2.5)
line.SetLineStyle(2)
line.SetLineColor(ROOT.kRed)
line.Draw("same,r")
can.Update()

#can.SaveAs("gif/limitVsLambdaS.gif")


outfile = TFile("out/limitPlot_%s.root"%theDir,"recreate")
outfile.cd()
limitVsMass.Write("limitVsMass")
outfile.Write()
outfile.Close()


#gres[0].Draw("ape")


'''


