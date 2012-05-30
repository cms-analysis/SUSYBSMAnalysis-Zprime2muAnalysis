#!/usr/bin/env python
import ROOT, os 
from ROOT import *

#from python.ufPlotTools import *
#setTDRStyle()

sbasedir = "/Users/kypreos/physics/zprime2012/code/EfficiencyResolutionFromMC/root3/" 
sbasedir = "root/" 

can = TCanvas("can","can",600,600)
ffit = TF1("fit","[0]*TMath::Voigt(x-[1],0,[2])")

gr = TGraph()
os.system('mkdir -p plots/zprimeGammaFits/')
for m in range(750,2301,250):
	
	infile= TFile(sbasedir+"crab_ana_effres_zp%d.root"%(m),"open")
	h1= infile.Get("EfficiencyFromMC/NumTotalTrigEff");

	mu,sigma =h1.GetMean(),h1.GetRMS() 
	ffit.SetParameters(h1.GetMaximum(),h1.GetMean(),h1.GetRMS())
#	ffit.FixParameter(2,1e-3)
	ffit.SetRange(m-2*sigma,mu+2*sigma)

	h1.Fit("fit","b")
	h1.Fit("fit","rb")
	ipoint = gr.GetN()
	gr.SetPoint(ipoint,m,ffit.GetParameter(2))
	
	print h1.GetRMS()
	
gr.SetTitle(";Z'_{#psi} mass [GeV];#Gamma [GeV]")
gr.Draw("ap")
gr.GetXaxis().SetNdivisions(105)

gr.Fit("pol1")
can.Update()
can.SaveAs("plots/zprimeGammaFits/zprimeWidthFits.png")
#can.Update()
