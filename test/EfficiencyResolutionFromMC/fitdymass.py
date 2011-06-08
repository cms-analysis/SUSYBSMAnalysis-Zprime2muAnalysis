#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.TH1.AddDirectory(0)

ps = plot_saver('plots/fitdymass')

int_lumi = 1000.
rebin = 2

masses = [20, 120, 200, 500, 800, 1000]
nevents = [2148325, 54550] + [55000] * 4
sigmas = [1631, 7.9, 0.97, 0.027, 0.0031, 9.7e-4]
for i in xrange(len(sigmas) - 1):
    sigmas[i] - sigmas[i+1]
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]

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
    ps.save('rawmass%i')
    hists.append(h)

htot = hists[0].Clone('htot')
htot.Sumw2()
for j in xrange(1, len(hists)):
    htot.Add(hists[j])
htot.SetTitle('')
htot.Draw()
htot.GetXaxis().SetTitle('reconstructed M _{#mu#mu} (GeV)')
htot.GetYaxis().SetTitle('Events/%i GeV' % rebin) # assumes original started out with 1 GeV bins
#htot.Scale(1/htot.Integral(htot.FindBin(300), htot.FindBin(2000)))
    
def fit_just_tail(lo, hi):
    fcn = ROOT.TF1('fcn', 'exp([0] + [1] * x**[2])', lo, hi)
    fcn.SetParLimits(0, 0, 10000)
    fcn.SetParLimits(1, -10, 0)
    fcn.SetParLimits(2, 0, 1)
    #fcn.FixParameter(2, 0.3)
    fcn.SetParNames("Norm", "a", "b")
    fcn.SetLineColor(ROOT.kBlue)
        
    htot.Fit(fcn, 'LVR')

    ps.c.Update()
    s = htot.GetListOfFunctions().FindObject("stats")
    s.SetOptFit(111)
    s.Draw()

    ps.save('mass%i_%i' % (lo, hi))

#l = [120] + range(200, 2000, 200) + [2000]
#for lo,hi in zip(l,l[1:]):
#    fit_just_tail(lo,hi)
#    fit_just_tail(lo, 2000)

fit_just_tail(300, 2000)
