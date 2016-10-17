#!/usr/bin/env python
#
#- Combines various sources of syst. uncertainty in the background shape
#- and parameterizes the uncertainty as a function of the invariant mass.
#

import os, sys
from array import array
from math import sqrt
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
#import fitdymass
from fitdymass import low,high,fitlow,fithigh
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetOptFit(0)

lg = ROOT.TLegend(0.15, 0.15, 0.5, 0.4)

ps = plot_saver('plots/bckgshape_syst')

mass_low  =   90.
mass_high = 3500.

#
#- Efficiency variation vs mass, from fit to simulation
#
eff_unc = ROOT.TF1('pol1', '[0] + x*[1]', mass_low, mass_high)
eff_unc.SetParameters(1.001, -1.2275e-5)
eff_unc.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
eff_unc.SetTitle('')
eff_unc.SetMinimum(0.50)
eff_unc.SetMaximum(1.01)
#eff_unc.SetMaximum(1.50)
lg.AddEntry(eff_unc,'Efficiency', 'LPE')
eff_unc.Draw()
#
#- Residual mass-dependent QCD NNLO/NLO corrections vs mass: assume 3% flat
#
qcd_corr = ROOT.TF1('pol0', '[0]', mass_low, mass_high)
qcd_corr.SetParameter(0, 0.97)
qcd_corr.SetLineColor(ROOT.kBlack)
qcd_corr.SetLineStyle(4)
lg.AddEntry(qcd_corr,'QCD Corr.','LPE')
qcd_corr.Draw("same")
#
#- EWK NLO corrections vs mass, from AN-2012/496
#
ewk_corr = ROOT.TF1('pol1', '[0] + x*[1]', mass_low, mass_high)
ewk_corr.SetParameters(1.01, -4.2e-5)
ewk_corr.SetLineColor(ROOT.kRed)
ewk_corr.SetLineStyle(2)
lg.AddEntry(ewk_corr,'EWK NLO Corr.','LPE')
ewk_corr.Draw("same")
#
#- Alignment uncertainty vs mass, from C1-C2 comparison
#
align_unc = ROOT.TF1('poly2', '[0] + x*[1] + x*x*[2]', mass_low, mass_high)
align_unc.SetParameters(0.9958, 1.259e-5, -1.453e-8)
align_unc.SetLineColor(ROOT.kBlue)
align_unc.SetLineStyle(4)
lg.AddEntry(align_unc,'Alignment Unc.','LPE')
align_unc.Draw("same")
#
#- Mass scale uncertainty vs mass: assume 5% per TeV
#
scale_corr = ROOT.TF1('pol1', '[0] + x*[1]', mass_low, mass_high)
scale_corr.SetParameters(1.01, -5.e-5)
scale_corr.SetLineColor(ROOT.kGray)
scale_corr.SetLineStyle(5)
lg.AddEntry(scale_corr,'Mass Scale Unc.','LPE')
scale_corr.Draw("same")
#
#- PDF uncertainties vs mass, e-mail from D. Bourilkov, 13/02/2013
#
masses = [  90,  250,  450,  950, 1450, 1950, 2450, 2950]
sigmas = [5.85, 4.64, 4.51, 7.06, 11.0, 15.6, 20.9, 31.4]
l = len(masses)
for i in range(1, l+1):
    sigmas[i-1] = (1. - sigmas[i-1]/100.)
print 'PDF uncertainties (2013 version):'
for i in range(1, l+1):
    print masses[i-1], sigmas[i-1]
#masses = array('d', masses)
#sigmas = array('d', sigmas)
#gpdf = ROOT.TGraph(l, masses, sigmas)
#pdf_unc = ROOT.TF1('pdf_unc', '[0] + x*[1] + x*x*[2]', mass_low, mass_high)
#gpdf.Fit(pdf_unc, "N")
#gpdf.SetLineColor(ROOT.kGreen)
#gpdf.SetMarkerColor(gpdf.GetLineColor())
#gpdf.SetMarkerStyle(20);
#gpdf.SetMarkerSize(1.2);
#gpdf.Draw("same P")
#pdf_unc.SetLineColor(gpdf.GetLineColor())
#pdf_unc.Draw("same")
#
#- PDF uncertainties vs mass, 2014 version, v2 of AN-12-348 (Fig. 3)
#
pdf_unc = ROOT.TF1('pdf_unc', '1. + [0] + x*[1] + x*x*[2]', mass_low, mass_high)
pdf_unc.SetParameters(-0.0276, -3.03e-5, -2.38e-8)
pdf_unc.SetLineColor(ROOT.kGreen)
lg.AddEntry(pdf_unc,'PDF Unc.','LPE')
pdf_unc.Draw("same")
masses = [  90,  250,  450,  950, 1450, 1950, 2450, 2950]
sigmas = [   0,    0,    0,    0,    0,    0,    0,    0]
l = len(masses)
print 'PDF uncertainties (2014 version):'
for i in range(1, l+1):
    print masses[i-1], pdf_unc.Eval(masses[i-1])
#
#- Total uncertainty
#
i = 0
mass = 100.
total_masses = []
total_sigmas = []
print '%10s%15s%15s%15s%15s%15s%15s%15s' % ('mass (GeV)', 'Efficiency', 'QCD NNLO', 'EWK NLO', 'Alignment', 'Momentum scale', 'PDF', 'Total')
#while mass < mass_high:
while mass < mass_high:
    total = sqrt((1.-eff_unc.Eval(mass))**2 + (1.-qcd_corr.Eval(mass))**2 + (1.-ewk_corr.Eval(mass))**2 + (1.-align_unc.Eval(mass))**2 + (1.-scale_corr.Eval(mass))**2 + (1.-pdf_unc.Eval(mass))**2)
    print '%10.f%15.4f%15.4f%15.4f%15.4f%15.4f%15.4f%15.4f' % (mass, 1.-eff_unc.Eval(mass), 1.-qcd_corr.Eval(mass), 1.-ewk_corr.Eval(mass), 1.-align_unc.Eval(mass), 1.-scale_corr.Eval(mass), 1.-pdf_unc.Eval(mass), total)
    total_masses.append(mass)
    total_sigmas.append(1.-total)
#    total_sigmas.append(1.+total)
    i = i + 1
    mass = mass + 100.
l = len(total_masses)
total_masses = array('d', total_masses)
total_sigmas = array('d', total_sigmas)
gtot = ROOT.TGraph(l, total_masses, total_sigmas)
gtot.SetMarkerStyle(21);
gtot.SetMarkerSize(1.2);
total_unc = ROOT.TF1('total_unc', '[0] + x*[1] + x*x*[2]', mass_low, mass_high)
gtot.Fit(total_unc, "N")
gtot.SetLineColor(ROOT.kMagenta)
gtot.SetMarkerColor(gtot.GetLineColor())
gtot.Draw("same P")
total_unc.SetLineColor(gtot.GetLineColor())
total_unc.Draw("same")
lg.AddEntry(total_unc,"Total Unc.",'LPE')
lg.Draw()
#
ps.save('unc_vs_mass', log=False)


#
# unc05 -- 5% at 3 TeV
# unc10 -- 10% at 3 TeV
# unc15 -- 15% at 3 TeV
# unc20 -- 20% at 3 TeV
# unc25 -- 25% at 3 TeV
# unc30 -- 30% at 3 TeV
# unc35 -- 35% at 3 TeV
# unc40 -- 40% at 3 TeV
# unc45 -- 45% at 3 TeV
# unc50 -- 50% at 3 TeV
#

lg2 = ROOT.TLegend(0.15,0.15,0.4,0.5)

unc05 = ROOT.TF1('unc05','[0] + [1]*x',low,fithigh)
unc05.SetParameters(1.0012,-1.71E-5)
unc05.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
#lg2.AddEntry(unc05,'5% at 3 TeV','L')

unc10 = ROOT.TF1('unc10','[0] + [1]*x',low,fithigh)
unc10.SetParameters(1.0024,-3.41E-5)
unc10.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
#lg2.AddEntry(unc10,'10% at 3 TeV','L')

unc15 = ROOT.TF1('unc15','[0] + [1]*x',low,fithigh)
unc15.SetParameters(1.0036,-5.12E-5)
unc15.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
lg2.AddEntry(unc15,'15% decrease at 3 TeV','L')

unc20 = ROOT.TF1('unc20','[0] + [1]*x',low,fithigh)
unc20.SetParameters(1.0048,-6.83E-5)
unc20.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
lg2.AddEntry(unc20,'20% decrease at 3 TeV','L')

unc25 = ROOT.TF1('unc25','[0] + [1]*x',low,fithigh)
unc25.SetParameters(1.0060,-8.53E-5)
unc25.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
lg2.AddEntry(unc25,'25% decrease at 3 TeV','L')

unc30 = ROOT.TF1('unc30','[0] + [1]*x',low,fithigh)
unc30.SetParameters(1.0072,-10.24E-5)
unc30.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
#lg2.AddEntry(unc30,'30% at 3 TeV','L')

unc35 = ROOT.TF1('unc35','[0] + [1]*x',low,fithigh)
unc35.SetParameters(1.0084,-11.94E-5)
unc35.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
#lg2.AddEntry(unc35,'35% at 3 TeV','L')

unc40 = ROOT.TF1('unc40','[0] + [1]*x',low,fithigh)
unc40.SetParameters(1.0096,-13.65E-5)
unc40.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
#lg2.AddEntry(unc40,'40% at 3 TeV','L')

unc45 = ROOT.TF1('unc40','[0] + [1]*x',low,fithigh)
unc45.SetParameters(1.0107,-15.36E-5)
unc45.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
#lg2.AddEntry(unc45,'45% at 3 TeV','L')

unc50 = ROOT.TF1('unc50','[0] + [1]*x',low,fithigh)
unc50.SetParameters(1.0119,-17.06E-5)
unc50.SetLineColor(ROOT.kBlue)
unc50.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
lg2.AddEntry(unc50,'50% decrease at 3 TeV','L')

# Plots of plus shapes
unc15plus = ROOT.TF1('unc15plus','2-unc15',low,fithigh)
unc20plus = ROOT.TF1('unc20plus','2-unc20',low,fithigh)
unc25plus = ROOT.TF1('unc25plus','2-unc25',low,fithigh)
unc50plus = ROOT.TF1('unc50plus','2-unc50',low,fithigh)

unc50plus.SetLineColor(ROOT.kOrange)
unc25plus.SetLineColor(ROOT.kGreen+2)
unc20plus.SetLineColor(ROOT.kBlue)
unc15plus.SetLineColor(ROOT.kBlack)
unc50plus.SetLineStyle(2)
unc25plus.SetLineStyle(2)
unc20plus.SetLineStyle(2)
unc15plus.SetLineStyle(2)
lg2.AddEntry(unc15plus,'15% increase at 3 TeV','L')
lg2.AddEntry(unc20plus,'20% increase at 3 TeV','L')
lg2.AddEntry(unc25plus,'25% increase at 3 TeV','L')
lg2.AddEntry(unc50plus,'50% increase at 3 TeV','L')

unc50.SetLineColor(ROOT.kOrange)
unc25.SetLineColor(ROOT.kGreen+2)
unc20.SetLineColor(ROOT.kBlue)
unc15.SetLineColor(ROOT.kBlack)
unc50.SetTitle('')
unc50.GetXaxis().SetRangeUser(100,fithigh)
unc50.Draw()
unc50.SetMinimum(0)
unc50.SetMaximum(2)
unc25.Draw('same')
unc20.Draw('same')
unc15.Draw('same')
unc50plus.Draw('same')
unc25plus.Draw('same')
unc20plus.Draw('same')
unc15plus.Draw('same')
lg2.Draw()
#ps.save('linear_tests',log=False)
