#!/usr/bin/python

# import ROOT in batch mode
import sys,os
import argparse
import math
import pickle

from setTDRStyle import setTDRStyle

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
from ROOT import TLegend
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ptbins = [52, 72, 100, 152, 200, 300, 452, 800.]

ROOT.gROOT.LoadMacro("cruijff.C+")
ROOT.gROOT.LoadMacro("doubleCB.C+")

rebinFactor = 2
xLow = 75
xHigh = 105

def getBinRange(histo,mlow,mhigh,reg="BB"):
	ymin =  0
	ymax = -1
	nbins = histo.GetNbinsY()
	for bin in range(nbins):
		if mlow==histo.GetYaxis().GetBinLowEdge(bin): 
			ymin = bin
		if mhigh==histo.GetYaxis().GetBinLowEdge(bin):
			ymax = bin

	if mhigh==ptbins[len(ptbins)-1]: ymax=-1
	if mhigh==ptbins[len(ptbins)-2] and "BE" in reg: ymax=-1
	return ymin,ymax

def loadHistos(inputdata,inputMC,region,weights,trackType,mcIsData,dataIsMC):
	_fileDATA = ROOT.TFile(inputdata)
	_fileMC   = []
	if mcIsData:
		_fileMC.append(ROOT.TFile(inputMC))
	else: 
		for mc in inputMC:
			_fileMC.append(ROOT.TFile(mc))
	
	hdata = ROOT.TH2F()
	hmc   = ROOT.TH2F()
	
	hdata.SetDirectory(0)
	hmc  .SetDirectory(0)
	ROOT.TH1.AddDirectory(ROOT.kFALSE)
	
	reg = ""
	if   ("BB" in region or "barrel" in region):  reg = "_BB"
	elif ("BE" in region or "endcap" in region):  reg = "_BE"
	hdata = _fileDATA.Get("Our2016MuonsPlusMuonsMinus%sResolution/DileptonMass_2d_vsPt%s" %(trackType,reg)).Clone()
	if mcIsData:
		hmc = _fileMC[0].Get("Our2016MuonsPlusMuonsMinus%sResolution/DileptonMass_2d_vsPt%s" %(trackType,reg)).Clone()
	else:	
		for k,mc in enumerate(inputMC):
			tmp   = _fileMC[k].Get("Our2016MuonsPlusMuonsMinus%sResolution/DileptonMass_2d_vsPt%s" %(trackType,reg)).Clone()
		if k==0 and not weights: 
			hmc = tmp
		elif k==0 and weights:
			print "Weighting with %s " %(weights[k])
			tmp.Scale(weights[k])
			hmc = tmp
		elif not weights:
			hmc.Add(tmp)
		else: 
			print "Weighting with %s " %(weights[k])
			tmp.Scale(weights[k])
			hmc.Add(tmp)
		print hmc.GetEntries()
	_fileDATA.Close()
	for f in _fileMC:
		f.Close()

	if "BB" in reg: 
		data = [ROOT.TH1D() for x in range(len(ptbins)-1)]
		mc   = [ROOT.TH1D() for x in range(len(ptbins)-1)]
		ptda = [0 for x in range(len(ptbins)-1)]
		ptmc = [0 for x in range(len(ptbins)-1)]
	else:
		data = [ROOT.TH1D() for x in range(len(ptbins)-2)]
		mc   = [ROOT.TH1D() for x in range(len(ptbins)-2)]
		ptda = [0 for x in range(len(ptbins)-2)]
		ptmc = [0 for x in range(len(ptbins)-2)]

	for h in data:
		h.SetDirectory(0)
		ROOT.TH1.AddDirectory(ROOT.kFALSE)
	for h in mc:
		h.SetDirectory(0)
		ROOT.TH1.AddDirectory(ROOT.kFALSE)
	
	for i,h in enumerate(data):        
		ymin,ymax=getBinRange(hdata,ptbins[i],ptbins[i+1],reg)
		hdata.GetYaxis().SetRangeUser(ptbins[i],ptbins[i+1])
		hmc.GetYaxis().SetRangeUser(ptbins[i],ptbins[i+1])
		ptda[i] = hdata.GetMean(2)
		ptmc[i] = hmc.GetMean(2)
		hdata.GetYaxis().SetRange()
		hmc  .GetYaxis().SetRange()

		data[i] = hdata.ProjectionX("datapy%s%s" %(ptbins[i],region),ymin,ymax)
		mc  [i] = hmc  .ProjectionX("mcpy%s%s" %(ptbins[i],region)  ,ymin,ymax)
		
		data[i].Rebin(rebinFactor)
		mc  [i].Rebin(rebinFactor)
		#~ data[i].Rebin(1)
		#~ mc  [i].Rebin(1)
		
		if (data[i].Integral() < 1500): 
			data[i].Rebin(4)
		if (mc[i].Integral() < 1500): 
			mc[i].Rebin(4)
		#~ if mcIsData:
			#~ mc[i].Rebin(2)
		#~ else:	
			#~ mc[i].Rebin(2)
		
#        if (ptbins[i]==200 or ptbins[i]==152) and "BE" in region:
#            mc[i].Rebin(2)

	return data,mc,ptda,ptmc

def doFit(hist,output,rap="BB",flavour="DATA",trackType="TunePNew"):
	c1 = ROOT.TCanvas("c1","c1",700,700)
	c1.cd()

	sig    = []
	sige   = []
	mean   = []
	meane  = []
	nChi2  = []
	
	DOCRYSTALBALL = False
	DOCRUIJFF = True
	DODOUBLECB = False
	
	for i,h in enumerate(hist):
		print "+++++++++++++++++++++++++++++++++++++++++"
		print "Fitting histogram for %d < pt_{l} <%d" %(ptbins[i],ptbins[i+1])
		print "+++++++++++++++++++++++++++++++++++++++++\n"

		fit_min = xLow
		fit_max = xHigh

		# fit with a gaussian 
		gaus = ROOT.TF1("gaus","gaus",fit_min,fit_max)
		gaus.SetLineColor(ROOT.kBlue)
		gaus.SetParameters(0,h.GetMean(),h.GetRMS())
		h.Fit("gaus","M0R+")
		print h.GetEntries()
		if DOCRYSTALBALL:
			funct = ROOT.TF1("crystal","crystalball",fit_min,fit_max)
			funct.SetLineColor(ROOT.kRed)
	   		ws = ROOT.RooWorkspace("tempWS")

			mass = ROOT.RooRealVar('mass','mass',91, fit_min, fit_max )
			getattr(ws,'import')(mass,ROOT.RooCmdArg())
			if h.Integral() < 1500:
				nDOF = (fit_max-fit_min)*2/(rebinFactor*4)-3
			else:	
				nDOF = (fit_max-fit_min)*2/rebinFactor-3
			dataHist = ROOT.RooDataHist("hist","hist",ROOT.RooArgList(ws.var("mass")),h)
			getattr(ws,'import')(dataHist,ROOT.RooCmdArg())
	#	    ws.factory("RooCBShape::cb(mass, mean[91.,80,110], sigma[2,0,10], alphaL[1,0,5], nL[1,0,5])")
			#~ ws.factory("RooCBShape::cb(mass, mean[0.0,-3,3], sigma[2,0,10], alphaL[3,-25,25], nL[5,-25,25])")
			ws.factory("RooCBShape::cb(mass, mean[0.0], sigma[2,0,10], alphaL[3,-25,25], nL[5,-25,25])")
			ws.factory("BreitWigner::bw(mass,meanZ[91.187], width[2.495])")
			bw = ws.pdf("bw")
			cb = ws.pdf("cb")
			mass.setBins(2000,"cache")
			mass.setMin("cache",0)
			mass.setMax("cache",1000); ## need to be adjusted to be higher than limit setting

			sigpdf = ROOT.RooFFTConvPdf("sig","sig",mass,bw,cb)
			getattr(ws,'import')(sigpdf,ROOT.RooCmdArg())

			fitResult = ws.pdf("sig").fitTo(ws.data("hist"),ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Minos(ROOT.kFALSE))
			#fitResult = ws.pdf("cb").fitTo(ws.data("hist"),ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Minos(ROOT.kFALSE))
			chi2 = ROOT.RooChi2Var("bla","blubb",ws.pdf("sig"),ws.data("hist")).getVal()
			print chi2
			nChi2.append(chi2/nDOF)
			#tmp_mean = gaus.GetParameter(1)
			#tmp_sig  = gaus.GetParameter(2)
			#min_mean = tmp_mean+tmp_mean
			#max_mean = tmp_mean-tmp_mean
			#funct.SetParameters(gaus.GetParameter(0), tmp_mean, tmp_sig, 1.4, 2)
			#funct.SetParLimits(0, 0, 2*gaus.GetParameter(0)) # //const
			#funct.SetParLimits(1, min_mean, max_mean)
			#funct.SetParLimits(2, 0, 2.5*h.GetRMS())# //sigma
			#funct.SetParLimits(3, 0.5, 2.)# //alpha
			#funct.SetParLimits(4, 0., 3.)# //alpha            
			#h.Fit("crystal","M0R+")
			funct.SetParameter(1,ws.var("meanZ").getVal())
			funct.SetParError(1,ws.var("meanZ").getError())
			funct.SetParameter(2,ws.var("sigma").getVal())
			funct.SetParError(2,ws.var("sigma").getError())
			funct.SetParameter(3,ws.var("alphaL").getVal())
			funct.SetParError(3,ws.var("alphaL").getError())
			funct.SetParameter(4,ws.var("nL").getVal())
			funct.SetParError(4,ws.var("nL").getError())
			# crystal
			mean  .append(ws.var("meanZ").getVal())
			meane .append(ws.var("meanZ").getError())
			sig   .append(ws.var("sigma").getVal())
			sige  .append(ws.var("sigma").getError())
		
		elif DOCRUIJFF:


			ROOT.gSystem.Load("./RooCruijff_cxx.so")
			ws = ROOT.RooWorkspace("tempWS")

			mass = ROOT.RooRealVar('mass','mass',91, fit_min, fit_max )
			getattr(ws,'import')(mass,ROOT.RooCmdArg())

			dataHist = ROOT.RooDataHist("hist","hist",ROOT.RooArgList(ws.var("mass")),h)
			getattr(ws,'import')(dataHist,ROOT.RooCmdArg())
	#	    ws.factory("RooCBShape::cb(mass, mean[91.,80,110], sigma[2,0,10], alphaL[1,0,5], nL[1,0,5])")
			#~ ws.factory("RooCruijff::cb(mass, mean[0.0,-3,3], sigma[2,0,10], sigma, alphaL[1,0,25], alphaR[1,0,25])")
			ws.factory("RooCruijff::cb(mass, mean[0.0], sigma[2,0,20], sigma, alphaL[1,0,25], alphaR[1,0,25])")

			if h.Integral() < 1500:
				nDOF = (fit_max-fit_min)*2/(rebinFactor*4)-3
			else:	
				nDOF = (fit_max-fit_min)*2/rebinFactor-3

			ws.factory("BreitWigner::bw(mass,meanZ[91.187], width[2.495])")
			bw = ws.pdf("bw")
			cb = ws.pdf("cb")
			mass.setBins(2000,"cache")
			mass.setMin("cache",0)
			mass.setMax("cache",1000); ## need to be adjusted to be higher than limit setting

			sigpdf = ROOT.RooFFTConvPdf("sig","sig",mass,bw,cb)
			getattr(ws,'import')(sigpdf,ROOT.RooCmdArg())

			fitResult = ws.pdf("sig").fitTo(ws.data("hist"),ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Minos(ROOT.kFALSE))
			chi2 = ROOT.RooChi2Var("bla","blubb",ws.pdf("sig"),ws.data("hist")).getVal()
			print chi2
			funct = ROOT.TF1("cruijff",ROOT.cruijff,fit_min,fit_max,5)
			funct.SetParNames("Constant","Mean","Sigma","AlfaL","AlfaR")
			#funct.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1),gaus.GetParameter(2),1.4,1.3)
			#funct.SetLineColor(ROOT.kRed)
			#h.Fit("cruijff","M0R+")
			funct.SetParameter(1,ws.var("meanZ").getVal())
			funct.SetParError(1,ws.var("meanZ").getError())
			funct.SetParameter(2,ws.var("sigma").getVal())
			funct.SetParError(2,ws.var("sigma").getError())
			funct.SetParameter(3,ws.var("alphaL").getVal())
			funct.SetParError(3,ws.var("alphaL").getError())
			funct.SetParameter(4,ws.var("alphaR").getVal())
			funct.SetParError(4,ws.var("alphaR").getError())
		   	nChi2.append(chi2/nDOF)
			mean  .append(funct.GetParameter(1))
			meane .append(funct.GetParError(1))
			sig   .append(funct.GetParameter(2))
			sige  .append(funct.GetParError(2))
		elif DODOUBLECB:


			ROOT.gSystem.Load("./RooDCBShape_cxx.so")

			ws = ROOT.RooWorkspace("tempWS")

			mass = ROOT.RooRealVar('mass','mass',91, fit_min, fit_max )
			getattr(ws,'import')(mass,ROOT.RooCmdArg())

			dataHist = ROOT.RooDataHist("hist","hist",ROOT.RooArgList(ws.var("mass")),h)
			getattr(ws,'import')(dataHist,ROOT.RooCmdArg())
			#~ ws.factory("RooDCBShape::cb(mass, mean[0.0,-3,3], sigma[2,0,20], alphaL[1,0,25] , alphaR[1,0,25], nL[5,0,25], nR[5,0,25])")
			ws.factory("RooDCBShape::cb(mass, mean[0.0], sigma[2,0,20], alphaL[1,0,25] , alphaR[1,0,25], nL[5,0,25], nR[5,0,25])")
			ws.factory("BreitWigner::bw(mass,meanZ[91.187], width[2.495])")
			bw = ws.pdf("bw")
			cb = ws.pdf("cb")
			mass.setBins(2000,"cache")
			mass.setMin("cache",0)
			mass.setMax("cache",1000); ## need to be adjusted to be higher than limit setting

			sigpdf = ROOT.RooFFTConvPdf("sig","sig",mass,bw,cb)
			getattr(ws,'import')(sigpdf,ROOT.RooCmdArg())

			fitResult = ws.pdf("sig").fitTo(ws.data("hist"),ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Minos(ROOT.kFALSE))

			chi2 = ROOT.RooChi2Var("bla","blubb",ws.pdf("sig"),ws.data("hist")).getVal()


			if h.Integral() < 1500:
				nDOF = (fit_max-fit_min)*2/(rebinFactor*4)-5
			else:	
				nDOF = (fit_max-fit_min)*2/rebinFactor-5


			funct = ROOT.TF1("doubleCB",ROOT.doubleCB,fit_min,fit_max,7)
			funct.SetParNames("Constant","Mean","Sigma","AlfaL","AlfaR","nL","nR")
			#funct.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1),gaus.GetParameter(2),ws.var("alphaL").getVal(),ws.var("alphaR").getVal(),ws.var("nL").getVal(),ws.var("nR").getVal())
			funct.SetParameters(gaus.GetParameter(0), ws.var("mean").getVal(),ws.var("sigma").getVal(),ws.var("alphaL").getVal(),ws.var("alphaR").getVal(),ws.var("nL").getVal(),ws.var("nR").getVal())
			funct.SetParError(1,ws.var("mean").getError())
			funct.SetParError(2,ws.var("sigma").getError())
			funct.SetParError(3,ws.var("alphaL").getError())
			funct.SetParError(4,ws.var("alphaR").getError())
			funct.SetParError(5,ws.var("nL").getError())
			funct.SetParError(6,ws.var("nR").getError())
			#funct.SetParLimits(3,1,3)
			#funct.SetParLimits(4,1,3)
			funct.SetLineColor(ROOT.kRed)
			#h.Fit("doubleCB","M0RL+")
				
			mean  .append(funct.GetParameter(1))
			meane .append(funct.GetParError(1))
			sig   .append(funct.GetParameter(2))
			sige  .append(funct.GetParError(2))
			nChi2.append(chi2/nDOF)
	
		   
		else:
			mean  .append(gaus.GetParameter(1))
			meane .append(gaus.GetParError(1))
			sig   .append(gaus.GetParameter(2))
			sige  .append(gaus.GetParError(2))
  
	
		plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
		style = setTDRStyle()
		ROOT.gStyle.SetOptStat(0)
		plotPad.UseCurrentStyle()
		plotPad.Draw()	
		plotPad.cd()
		ROOT.gStyle.SetTitleXOffset(1.45)	

		if DODOUBLECB or DOCRYSTALBALL or DOCRUIJFF:
			ws.var("mass").setBins(30)
			frame = ws.var('mass').frame(ROOT.RooFit.Title('Invariant mass of dimuon pairs'))
			frame.GetXaxis().SetTitle('m_{#mu#mu} [GeV]')
			frame.GetYaxis().SetTitle("Events / 2 GeV")
			ROOT.RooAbsData.plotOn(ws.data('hist'), frame,ROOT.RooFit.Name("hist"))
			ws.pdf('sig').plotOn(frame,ROOT.RooFit.Name("sig"))
			frame.Draw()
			#chi2 = frame.chiSquare("sig","hist",nDOF) 
		else:

			h.GetXaxis().SetTitle("m_{ll} [GeV]")
			h.SetLineColor(ROOT.kBlack)
			h.GetXaxis().SetRangeUser(fit_min,fit_max)
			h.SetMarkerStyle(20)
			h.SetMarkerSize(0.7)
				
			h.Draw("E")
			if DOCRYSTALBALL or DOCRUIJFF or DODOUBLECB:
				funct.Draw("SAME")
			else:
				gaus.Draw("SAME")
			
		latex = ROOT.TLatex()
		latex.SetTextFont(42)
		latex.SetTextAlign(31)
		latex.SetTextSize(0.04)
		latex.SetNDC(True)
		latexCMS = ROOT.TLatex()
		latexCMS.SetTextFont(61)
		latexCMS.SetTextSize(0.055)
		latexCMS.SetNDC(True)
		latexCMSExtra = ROOT.TLatex()
		latexCMSExtra.SetTextFont(52)
		latexCMSExtra.SetTextSize(0.03)
		latexCMSExtra.SetNDC(True)

		latex.DrawLatex(0.95, 0.96, "(13 TeV)")

		cmsExtra = "Preliminary" 
		latexCMS.DrawLatex(0.78,0.88,"CMS")
		yLabelPos = 0.84
		latexCMSExtra.DrawLatex(0.78,yLabelPos,"%s"%(cmsExtra))

		latexFit1 = ROOT.TLatex()
		latexFit1.SetTextFont(61)
		latexFit1.SetTextSize(0.035)
		latexFit1.SetNDC(True)
		latexFit1.DrawLatex(0.19, 0.84, "%d < p_{T} <%d" %(ptbins[i],ptbins[i+1]))
		
		latexFit = ROOT.TLatex()
		latexFit.SetTextFont(42)
		latexFit.SetTextSize(0.030)
		latexFit.SetNDC(True)        
		latexFit.DrawLatex(0.19, 0.8,"%s = %5.3g #pm %5.3g"%("mean bias",ws.var("mean").getVal(),ws.var("mean").getError()))
		for par in range(funct.GetNpar()-1):
			yPos = 0.74-0.04*(float(par))
			latexFit.DrawLatex(0.19, yPos,"%s = %5.3g #pm %5.3g"%(funct.GetParName(par+1),funct.GetParameter(par+1),funct.GetParError(par+1)))
		if funct.GetNDF()> 0:
				latexFit.DrawLatex(0.19, 0.54, "#chi^{2}/ndf = %5.1f / %2.0f = %4.2f" %(funct.GetChisquare(),funct.GetNDF(),funct.GetChisquare()/funct.GetNDF()))

		if DODOUBLECB or DOCRYSTALBALL or DOCRUIJFF:
				latexFit.DrawLatex(0.19, 0.54, "#chi^{2}/ndf = %5.1f / %2.0f = %4.2f" %(chi2,nDOF,chi2/nDOF))
		saveas = "/MassRes_%s_%s_Pt%d_%d_%s" %(trackType,flavour,ptbins[i],ptbins[i+1],rap)
		c1.SaveAs(output+saveas+".root")
		c1.SaveAs(output+saveas+".C")
		c1.SaveAs(output+saveas+".png")
		c1.SaveAs(output+saveas+".pdf")

	print "DONE Fitting..."
	return mean,meane,sig,sige, nChi2

	
def drawMassRes(data,mc,output,rapidity,ptda,ptmc,trackType,mcIsData,dataIsMC):
	style = setTDRStyle()
	
	pt_e = [0 for x in range(len(data))]
	pt_x = [0 for x in range(len(data))]
	for i,pt in enumerate(pt_x):
		pt_x[i] = ptbins[i]+(ptbins[i+1]-ptbins[i])/2.
		pt_e[i] = (ptbins[i+1]-ptbins[i])/2.
	if dataIsMC:
		(da_mean,da_meane,da_sig,da_sige, da_nChi2) = doFit(data,output,rapidity,"MC2",trackType)
	else:	
		(da_mean,da_meane,da_sig,da_sige, da_nChi2) = doFit(data,output,rapidity,"DATA",trackType)
	if mcIsData:	
		(mc_mean,mc_meane,mc_sig,mc_sige, mc_nChi2) = doFit(mc  ,output,rapidity,"DATA2",trackType)
	else:
		(mc_mean,mc_meane,mc_sig,mc_sige, mc_nChi2) = doFit(mc  ,output,rapidity,"MC",trackType)
	result = {}
	result["data_mean"] = da_mean
	result["data_meane"] = da_meane
	result["data_sig"] = da_sig
	result["data_sige"] = da_sige
	result["mc_mean"] = mc_mean
	result["mc_meane"] = mc_meane
	result["mc_sig"] = mc_sig
	result["mc_sige"] = mc_sige
	result["ptda"] = ptda
	result["ptmc"] = ptmc
	result["da_nChi2"] = da_nChi2
	result["mc_nChi2"] = mc_nChi2


	print da_sig, mc_sig

	pklFile = open(output+"/MassResolutionVsPt_%s_%s.pkl" %(trackType,rapidity),"w")
	pickle.dump(result,pklFile)
	pklFile.close()
	
	c2 = ROOT.TCanvas("c2","c2",700,700)
	c2.cd()

	# Upper plot will be in pad1
	pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetGrid()        # Vertical grid
	pad1.SetBottomMargin(0)
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad
	pad1.SetTicks()
  
	res_data  = ROOT.TGraphAsymmErrors(len(pt_x))
	res_data.SetName("res_data")
	res_mc    = ROOT.TGraphAsymmErrors(len(pt_x))
	res_mc  .SetName("res_mc")
	ratio     = ROOT.TGraphErrors(len(pt_x))
	ratio   .SetName("ratio")
	print len(pt_x)
	
	for i,pt in enumerate(pt_x):        
		res_data.SetPoint(i,ptda[i],da_sig[i])
		res_data.SetPointError(i,ptda[i]-ptbins[i],ptbins[i+1]-ptda[i],da_sige[i],da_sige[i])
		res_mc  .SetPoint(i,ptmc[i],mc_sig[i])
		res_mc  .SetPointError(i,ptmc[i]-ptbins[i],ptbins[i+1]-ptmc[i],mc_sige[i],mc_sige[i])
		if mc_sig[i] > 0:
			ratio   .SetPoint(i,pt,da_sig[i]/mc_sig[i])
			ratio   .SetPointError(i,pt_e[i],(da_sig[i]/mc_sig[i])*math.sqrt((da_sige[i]/da_sig[i])**2+(mc_sige[i]/mc_sig[i])**2))
	res_data.SetMarkerStyle(22)
	res_data.SetMarkerColor(ROOT.kBlack)
	res_data.SetLineColor(ROOT.kBlack)
	res_data.SetFillColor(0)
	res_data.SetTitle("Dimuon mass resolution vs pT for %s tracks"%trackType)
	res_data.GetYaxis().SetTitle("#sigma (Z peak)")
	res_data.GetYaxis().SetTitleOffset(1.2)
	res_data.GetYaxis().SetRangeUser(0.,6.)
	if trackType == "Outer":
		res_data.GetYaxis().SetRangeUser(1.,20.)
	res_data.GetXaxis().SetRangeUser(ptbins[0],ptbins[len(ptda)])
	res_data.Draw("AP E0")
		
	res_mc.SetMarkerStyle(22)
	res_mc.SetMarkerColor(ROOT.kRed)
	res_mc.SetLineColor(ROOT.kRed)
	res_mc.SetFillColor(0)
	res_mc.SetTitle("Dimuon mass resolution vs pT for %s tracks"%trackType)
	res_mc.GetYaxis().SetTitle("Resolution (Z peak)")
	res_mc.GetYaxis().SetTitleOffset(1.5)
	res_mc.Draw("P E0 SAME")
	leg = TLegend(.35,.7,.50,.80,"","brNDC")
	if mcIsData:
		leg.AddEntry(res_data,"DATA 2017")
		leg.AddEntry(res_mc,"DATA 2016")
	elif dataIsMC:
		leg.AddEntry(res_data,"MC 2017")
		leg.AddEntry(res_mc,"MC 2016")		
	else:
		leg.AddEntry(res_data,"DATA","p")
		leg.AddEntry(res_mc,"Simulation")
	leg.SetTextFont(42)
	leg.SetBorderSize(0)
	leg.SetTextSize(.04)
	leg.Draw("SAME")
	latex = ROOT.TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = ROOT.TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.055/0.7)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.03/0.7)
	latexCMSExtra.SetNDC(True)
	
	latex.DrawLatex(0.95, 0.96, "(13 TeV)")
	
	cmsExtra = "Preliminary" 
	latexCMS.DrawLatex(0.78,0.88,"CMS")
	yLabelPos = 0.84
	latexCMSExtra.DrawLatex(0.78,yLabelPos,"%s"%(cmsExtra))
	c2.cd()          # Go back to the main canvas before defining pad2
	pad2 = ROOT.TPad("pad2", "pad2",0, 0.1, 1, 0.30)    
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.25)
	pad2.SetGrid()
	pad2.Draw()
	pad2.cd()
	pad2.SetTicks()

	ratio.SetMarkerColor(ROOT.kBlue-4)
	ratio.SetFillColor(ROOT.kBlue-4 )
	ratio.SetTitle("")
	ratio.GetYaxis().SetTitle("Data/MC")
	if mcIsData:
		ratio.GetYaxis().SetTitle("Data 2017 / Data 2016")
	elif dataIsMC:
		ratio.GetYaxis().SetTitle("MC 2017 / MC 2016")
	ratio.GetXaxis().SetTitle("p_{T} (#mu^{#pm}) [GeV]")
	ratio.GetYaxis().SetRangeUser(0.5,1.5)
	ratio.GetXaxis().SetRangeUser(ptbins[0],ptbins[len(ptda)])
	ratio.GetYaxis().SetTitleOffset(0.50)
	ratio.GetYaxis().SetTitleSize(0.14)
	ratio.GetYaxis().SetLabelSize(0.14)    
	ratio.GetYaxis().SetNdivisions(506)    
	ratio.GetXaxis().SetTitleSize(0.2)
	ratio.GetXaxis().SetTitleOffset(1.2)
	ratio.GetXaxis().SetLabelSize(0.20)
			
	ratio.Draw("A P E2")
	pad2.Update()
	line = ROOT.TLine(ptbins[0],1,ptbins[len(ptda)],1)

	line.SetLineColor(ROOT.kBlue+1)
	line.SetLineWidth(2)
	line.Draw()
	saveas = "/MassResolutionVsPt_%s_%s" %(trackType,rapidity)
	c2.SaveAs(output+saveas+".png")
	c2.SaveAs(output+saveas+".pdf")
	c2.SaveAs(output+saveas+".root")
	c2.SaveAs(output+saveas+".C")
	
	# PRINT FIT RESULTS!!!
	print "|--------------------------------------------------------------------------|"
	print "|                   MASS RESOLUTION PARAMETRIZATION                        |"
	print "|--------------------------------------------------------------------------|"
	print "|      pT   %s      |      Mean:  Data / MC      |    Sigma: Data / MC     | " %(rapidity)
	print "|--------------------------------------------------------------------------|"
	
	for i in range(0,len(ptda)):
		if ptbins[i]>200 and rapidity=="E": break
	if mc_mean[i] > 0:
		print "| %4d < p_T < %4d |  %5.3f / %5.3f ( %5.3f ) | %5.3f / %5.3f ( %5.3f ) |" %(pt_x[i]-pt_e[i], pt_x[i]+pt_e[i],
																					   da_mean[i], mc_mean[i], da_mean[i]/mc_mean[i],
																						   da_sig[i], mc_sig[i], da_sig[i]/mc_sig[i])
	print "|--------------------------------------------------------------------------|"

	
def makeMassRes(inputDATA,inputMC,output,weights,trackType,mcIsData,dataIsMC):
	style = setTDRStyle()
	ROOT.gStyle.SetTitleYOffset(1.45)
	ROOT.gStyle.SetTitleXOffset(1.45)
	ROOT.gStyle.SetOptFit(0)
	ROOT.gStyle.SetStatX(.9)
	ROOT.gStyle.SetStatY(.9)
	
	(data_B,mc_B,ptdaB,ptmcB) = loadHistos(inputDATA,inputMC,"BB",weights,trackType,mcIsData,dataIsMC)
	(data_E,mc_E,ptdaE,ptmcE) = loadHistos(inputDATA,inputMC,"BE",weights,trackType,mcIsData,dataIsMC)

	drawMassRes(data_B,mc_B,output,"BB",ptdaB,ptmcB,trackType,mcIsData,dataIsMC)
	drawMassRes(data_E,mc_E,output,"BE",ptdaE,ptmcE,trackType,mcIsData,dataIsMC)
	
#    
#    resBB.SetMarkerStyle(22)
#    resBB.SetMarkerColor(ROOT.kRed)
#    resBB.SetLineColor(ROOT.kRed)
#    resBB.SetFillColor(0)
#    resBB.SetTitle("Dimuon mass resolution vs pT")
#    resBB.GetYaxis().SetTitle("Dimuon Mass Resolution")
#    resBB.GetYaxis().SetTitleOffset(1.5)
#    resBB.GetXaxis().SetTitle("p_T (#mu^{#pm}) [GeV]")
#    resBB.GetYaxis().SetRangeUser(0,.2)
#    resBB.GetXaxis().SetRangeUser(0,5000)
#    resBB.GetFunction("fun").SetLineColor(ROOT.kRed+1)
#    resBB.Draw("AP E0")
#    
#    resBE.SetMarkerStyle(22)
#    resBE.SetMarkerColor(ROOT.kBlue+1)
#    resBE.SetLineColor(ROOT.kBlue+1)
#    resBE.SetFillColor(0)
#    resBE.SetTitle("Dimuon mass resolution vs mass")
#    resBE.GetYaxis().SetTitle("Dimuon Mass Resolution")
#    resBE.GetYaxis().SetTitleOffset(1.5)
#    resBE.GetXaxis().SetTitle("p_T (#mu^{#pm}) [GeV]")
#    resBE.GetYaxis().SetRangeUser(0,.2)
#    resBE.GetXaxis().SetRangeUser(0,5000)
#    resBE.GetFunction("fun").SetLineColor(ROOT.kBlue+2)
#    resBE.Draw("PE0 SAME")
#        
#    leg = ROOT.TLegend(.35,.7,.50,.80,"","brNDC")
#    leg.AddEntry(resBB,"BB")
#    leg.AddEntry(resBE,"BE")
#    leg.SetTextFont(42)
#    leg.SetBorderSize(0)
#    leg.SetTextSize(.02)
#    leg.Draw("SAME")
#    
#    res.SetGrid()
#    saveas = "/MassResolutionVsPt_2CAT"
#    res.SaveAs(output+saveas+".png")
#    res.SaveAs(output+saveas+".pdf")
#    res.SaveAs(output+saveas+".root")
#    res.SaveAs(output+saveas+".C")
	
		 
#### ========= MAIN =======================
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(usage="makeMassRes.py [options]",description="Compute mass resolution",
									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("--iDATA", dest="inputDATA",default="", help='Input filename')
	parser.add_argument("--iDATA2", dest="inputDATA2",default="", help='Input filename 2')
#    parser.add_argument("--iMC", dest="inputMC",default="", help='Input filename')
	parser.add_argument("--iMC", dest="inputMC", type=str, help='The list of input files, comma separated if more than one file',default = "")
	parser.add_argument("--iMC2", dest="inputMC2", type=str, help='The list of input files, comma separated if more than one file',default = "")
	parser.add_argument("--weight",type=str, help='The list of weights, comma separated if more than one',nargs=1)
	parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name to store results')
	args = parser.parse_args()
	
	if not os.path.exists(args.output):
		os.makedirs(args.output);
		if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+args.output)
	
	inputDATA = args.inputDATA
	inputDATA2 = args.inputDATA2
	mcIsData = False
	dataIsMC = False
	if args.inputMC == "":
		inputMC = inputDATA2
		mcIsData = True
	else:
		inputMC   = args.inputMC.split(",")
	if inputDATA == "":
		inputDATA = args.inputMC2
		dataIsMC = True
		
	print inputMC
	output=args.output

	weights = []     
	if args.weight is not None:
		tmp = args.weight[0].split(",")
		for i in tmp: 
			weights.append(float(i))
		print weights
		
	print "Running on: %s %s with 2 categories" %(inputDATA,inputMC)
	print "Saving result in: %s" %(output)

	tracks = ["Inner","Outer","Global","TPFMS","Picky","DYT","TunePNew"]
	#tracks = ["Outer","Global","TPFMS","Picky","DYT","TunePNew"]
	for trackType in tracks:	
		makeMassRes(inputDATA,inputMC,output,weights,trackType,mcIsData,dataIsMC)
	print "DONE"
