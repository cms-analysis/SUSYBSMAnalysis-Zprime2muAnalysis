#!/usr/bin/env py

import ROOT, commands, math
from ROOT import *

from ufPlotTools import *
setTDRStyle()

theDir = "out400b"

files = commands.getoutput("ls -l %s | awk {'print $9'}"%theDir).rsplit("\n")
files.remove(files[0]) #remove the empty first part

def DumpFiles(files):
	for ifile in files:
		print ifile

DumpFiles(files)


from ROOT import TFile,TGraphErrors,TF1

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





