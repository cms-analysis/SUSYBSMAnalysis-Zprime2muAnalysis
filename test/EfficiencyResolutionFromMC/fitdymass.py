#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.TH1.AddDirectory(0)

ps = plot_saver('plots/fitdymass')

int_lumi = 1000.
rebin = 5

masses = [20, 120, 200, 500, 800, 1000]
nevents = [2148325, 54550] + [55000] * 4
sigmas = [1631, 7.9, 0.97, 0.027, 0.0031, 9.7e-4] # in pb
for i in xrange(len(sigmas) - 1):
    sigmas[i] - sigmas[i+1]
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
#weights = [x/weights[-1] for x in weights]

hists = []
for m,w in zip(masses, weights):
    fn = 'ana_datamc_dy%i.root' % m if m != 20 else 'ana_datamc_zmumu.root'
    f = ROOT.TFile(fn)
    d = f.OurNewMuonsPlusMuonsMinusHistos
    
    h = d.Get('DileptonMass').Clone('dy%i' % m)
    h.Sumw2()
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(40, 2000)
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%i' % m)
    hists.append(h)

htot = hists[0].Clone('htot')
htot.Sumw2()
for j in xrange(1, len(hists)):
    htot.Add(hists[j])
htot.SetTitle('')
htot.Draw()
htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the in
#htot.Scale(1/htot.Integral(htot.FindBin(300), htot.FindBin(2000)))

#file = ROOT.TFile('dy.root', 'recreate')
#htot.SetDirectory(file)
#file.Write()
#file.Close
#raise 'done'
    
def fit_it(lo, hi):
 #   fcn = ROOT.TF1('fcn', 'exp([0] + [1] * x**[2])', lo, hi)
    fcn = ROOT.TF1('fcn', 'exp([0] + [1]*x)/x**[2]', lo, hi)
    fcn.SetParLimits(0, 0, 1000)
    fcn.SetParLimits(1, -20, 0)
    fcn.SetParLimits(2, 0, 10)
    fcn.SetParNames("N", "a", "b")
    fcn.SetLineColor(ROOT.kBlue)

    htot.Fit(fcn, 'LVR')

    ps.c.Update()
    s = htot.GetListOfFunctions().FindObject("stats")
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(10)
    s.SetOptFit(111)
    s.Draw()

    ps.save('mass%i_%i' % (lo, hi))

    xax = htot.GetXaxis()
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), 'relative residual;m(#mu^{+}#mu^{-}) (GeV);f(m_{i})/h_{i} - 1', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
    for h in [hres]:
        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= lo and xhi <= hi:
            res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
            if htot.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/htot.GetBinContent(i))
                hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))

    hres.Draw('e')
    ps.save('res_%i_%i' % (lo, hi), log=False)
    
                              
l = range(200, 2200, 200)
l = [400, 600, 800]
for lo in l:
    fit_it(lo, 2000)
#fit_it(400,2000)
