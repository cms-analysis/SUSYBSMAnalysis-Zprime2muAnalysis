#!/usr/bin/env python

import os
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)

os.system('mkdir -p plots/resfromzp')
c = ROOT.TCanvas('c', '', 820, 630)

masses = [1000, 1250, 1500, 1750]
sigma = []
esigma = []
chnam, val, err, xlolim, xhilim, iuint = ROOT.TString(''), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Long(1)

for m in masses:
    f = ROOT.TFile('ana_effres_zp%i.root' % m)
    d = f.respmc
    h = d.Get('DileptonMassRes')
    fcn = ROOT.TF1('mg', 'gaus', h.GetMean() - 1.5*h.GetRMS() , h.GetMean() + 1.5*h.GetRMS() )
    h.Fit(fcn, 'IVR')
    c.SaveAs('plots/resfromzp/mass_res_zp%i.png' % m)
    ROOT.gMinuit.mnpout(2, chnam, val, err, xlolim, xhilim, iuint)
    print val, err, chnam
    sigma.append(float(val))
    esigma.append(float(err))
    print

l = len(masses)
masses = array('d', masses)
emasses = array('d', l*[10.])
sigma = array('d', sigma)
esigma = array('d', esigma)

g = ROOT.TGraphErrors(l, masses, sigma, emasses, esigma)
g.Fit("pol2")
g.GetXaxis().SetTitle('dilepton mass (GeV)')
g.GetYaxis().SetTitle('dilepton relative mass resolution')
g.GetYaxis().SetLabelSize(0.02)
g.SetTitle('')
g.Draw("AP")
c.Update()
s = g.GetListOfFunctions().FindObject("stats")
s.SetOptFit(111)
s.SetY1NDC(0.15)
s.SetY2NDC(0.40)
s.SetX1NDC(0.50)
s.SetX2NDC(0.95)
c.Update()
c.SaveAs('plots/resfromzp/mass_res_zp.png')
c.SaveAs('plots/resfromzp/mass_res_zp.root')
