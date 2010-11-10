#!/usr/bin/env python

from math import exp
import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.TH1.AddDirectory(0)

c = ROOT.TCanvas('c', '', 820, 630)
c.SetLogy(1)
os.system('mkdir -p plots/fitdymass')

# 166.9 nb in jul15.root, 2880.2 in prompt.root
int_lumi = (166.9 + 2880.2)/1000 # /pb

masses = [20, 120, 200, 500, 800]
nevents = [1768457] + [55000] * 4
sigmas = [1686, 7.9, 0.965, 0.0269, 0.0031]
for i in xrange(len(sigmas) - 1):
    sigmas[i] - sigmas[i+1]
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]

hists = []
for m,w in zip(masses, weights):
    fn = 'ana_effres_dy%i.root' % m if m != 20 else 'ana_effres_zmumu.root'
    f = ROOT.TFile(fn)
    d = f.histospmc
    
    h = d.Get('DileptonMass').Clone('dy%i' % m)
    h.Rebin(3)
    h.GetXaxis().SetRangeUser(40, 400)
    h.Scale(w)
    h.Draw()
    c.SaveAs('plots/fitdymass/rawmass%i.png' % m)
    hists.append(h)

htot = hists[0].Clone('htot')
for j in xrange(1, len(hists)):
    htot.Add(hists[j])
htot.SetTitle('')
htot.Draw()
htot.GetXaxis().SetTitle('reconstructed M _{#mu#mu} (GeV)')
htot.GetYaxis().SetTitle('Events/10 GeV')

def fit_full_zlineshape():
    import z0dyshape
    fmix = ROOT.TF1("fmix", z0dyshape.MixFunc, 40, 400, 11);
    fmix.SetNpx(10000)
    fmix.FixParameter(0,  2.44298e+01);
    fmix.FixParameter(1,  6.99626e+00);
    fmix.FixParameter(2,  2.03121e-01);
    fmix.FixParameter(3,  9.45741e+06);
    fmix.FixParameter(4,  2.63045e+00);
    fmix.FixParameter(5,  1.80690e-02);
    fmix.FixParameter(6,  5.12014e+01);
    fmix.FixParameter(7,  0.00000e+00);
    fmix.FixParameter(8,  0.00000e+00);
    fmix.FixParameter(9,  2.19753e-02);
    fmix.SetParameter(10, 1);

    htot.SetStats(0)
    htot.Fit(fmix, 'ILVR')
    c.SaveAs('plots/fitdymass/full.root')
    c.SaveAs('plots/fitdymass/full.png')

def fit_just_tail(lo, hi):
    def f(x,p):
        #return p[0] * exp(p[1] * x[0]**0.3)
        return p[0] * exp(p[1] * x[0]**p[2])

    fcn = ROOT.TF1('fcn', f, lo, hi, 3)
    fcn.SetParLimits(0, 0, 1e10)
    fcn.SetParLimits(1, -10, 0)
    fcn.SetParLimits(2, 0, 1)
    fcn.SetParNames("Norm", "a", "b")
    #fcn.FixParameter(2, 0.3)

    htot.Fit(fcn, 'ILVR')

    c.Update()
    s = htot.GetListOfFunctions().FindObject("stats")
    s.SetOptFit(111)
    s.Draw()

    fn = 'plots/fitdymass/mass%i_%i.png' % (lo, hi)
    c.SaveAs(fn)
    c.SetLogy(1)
    c.SaveAs(fn.replace('.png', '.log.png'))
    c.SaveAs(fn.replace('.png', '.root'))

#l = [120] + range(200, 2000, 200) + [2000]
#for lo,hi in zip(l,l[1:]):
#    fit_just_tail(lo,hi)

fit_just_tail(300, 2000)
fit_full_zlineshape()
