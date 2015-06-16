#!/usr/bin/env python
#
#- Combines various sources of syst. uncertainty in the background shape
#- and parameterizes the uncertainty as a function of the invariant mass.

import os, sys
from array import array
from math import sqrt
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetOptFit(0)

ps = plot_saver('plots/bckgshape_syst')

mass_low  =   90.
mass_high = 3500.

#
#- Efficiency variation vs mass, from fit to simulation
#
eff_unc = ROOT.TF1('pol1', '[0] + x*[1]', mass_low, mass_high)
eff_unc.SetParameters(1.001, -1.2275e-5)
eff_unc.GetXaxis().SetTitle('dilepton mass (GeV)')
eff_unc.SetTitle('')
eff_unc.SetMinimum(0.50)
eff_unc.SetMaximum(1.01)
#eff_unc.SetMaximum(1.50)
eff_unc.Draw()
#
#- Residual mass-dependent QCD NNLO/NLO corrections vs mass: assume 3% flat
#
qcd_corr = ROOT.TF1('pol0', '[0]', mass_low, mass_high)
qcd_corr.SetParameter(0, 0.97)
qcd_corr.SetLineColor(ROOT.kBlack)
qcd_corr.SetLineStyle(4)
qcd_corr.Draw("same")
#
#- EWK NLO corrections vs mass, from AN-2012/496
#
ewk_corr = ROOT.TF1('pol1', '[0] + x*[1]', mass_low, mass_high)
ewk_corr.SetParameters(1.01, -4.2e-5)
ewk_corr.SetLineColor(ROOT.kRed)
ewk_corr.SetLineStyle(2)
ewk_corr.Draw("same")
#
#- Alignment uncertainty vs mass, from C1-C2 comparison
#
align_unc = ROOT.TF1('poly2', '[0] + x*[1] + x*x*[2]', mass_low, mass_high)
align_unc.SetParameters(0.9958, 1.259e-5, -1.453e-8)
align_unc.SetLineColor(ROOT.kBlue)
align_unc.SetLineStyle(4)
align_unc.Draw("same")
#
#- Mass scale uncertainty vs mass: assume 5% per TeV
#
scale_corr = ROOT.TF1('pol1', '[0] + x*[1]', mass_low, mass_high)
scale_corr.SetParameters(1.01, -5.e-5)
scale_corr.SetLineColor(ROOT.kGray)
scale_corr.SetLineStyle(5)
scale_corr.Draw("same")
#
#- PDF uncertainties vs mass, 2013 version, e-mail from D. Bourilkov, 13/02/2013
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
#
ps.save('unc_vs_mass', log=False)
