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

masses  = [     20,   120,   200,    500,     800,    1000,    1500,    2000]
nevents = [3293740, 99984, 99990,  99992,   99984,   99989,   99992,   99974]
#sigmas  = [  1915.,  12.2,  1.53, 0.0462, 0.00586, 0.00194, 1.70e-4, 2.21e-5] # in pb, PYTHIA*1.3
sigmas  = [  1915.,  12.2,  1.52, 0.0452, 0.00562, 0.00184, 1.75e-4, 2.26e-5] # in pb, POWHEG*1.024
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
#weights = [x/weights[-1] for x in weights]

hists = []
hists_dir = '../DataMCSpectraComparison/mc/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_dy%i_c1.root' % m if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir + fn
    f = ROOT.TFile(fn)
    d = f.Our2012MuonsPlusMuonsMinusHistos
    
    h = d.Get('DileptonMass').Clone('dy%i' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(40, 2500)
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%i' % m)
    hists.append(h)

htot = hists[0].Clone('htot')
for j in xrange(1, len(hists)):
    htot.Add(hists[j])
htot.SetTitle('')
htot.Draw()
htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
#htot.Scale(1/htot.Integral(htot.FindBin(300), htot.FindBin(2000)))

#file = ROOT.TFile('dy.root', 'recreate')
#htot.SetDirectory(file)
#file.Write()
#file.Close
#raise 'done'
    
def fit_it(lo, hi):
#    fcn = ROOT.TF1('fcn', 'exp([0] + [1] * x**[2])', lo, hi)
    fcn = ROOT.TF1('fcn', 'exp([0] + [1]*x)*x**[2]', lo, hi)
#    fcn = ROOT.TF1('fcn', 'exp([0] + [1]*x + [2]*x*x)*x**[3]', lo, hi)
#    fcn.SetParLimits(0, 0, 1000)
#    fcn.SetParLimits(1,  -1, 1)
#    fcn.SetParLimits(2, -10, 0)
    fcn.SetParNames("N", "a", "b")
#    fcn.SetParNames("N", "a", "b", "c")
    fcn.SetLineColor(ROOT.kBlue)

#    htot.Fit(fcn, 'LVR')
    htot.Fit(fcn, 'LVRW')

    ps.c.Update()
    s = htot.GetListOfFunctions().FindObject("stats")
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(10)
    s.SetOptFit(111)
    s.Draw()

    ps.save('mass%i_%i' % (lo, hi))

    xax = htot.GetXaxis()
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-hist)/hist', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
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

    hres.SetMinimum(-1.5)
    hres.SetMaximum( 1.5)
    hres.Draw('e')
    l1 = ROOT.TLine(lo, 0., hi,  0.)
    l1.Draw()
    ps.save('res_%i_%i' % (lo, hi), log=False)

def draw_overlay(fsets):
    ROOT.gStyle.SetOptStat(0);
    ROOT.gStyle.SetOptFit(1111);
    xmin = 100
    xmax = 2600
    rebin = 60
    htots = []
    for fset in fsets:
        hists = []
        for m,w in zip(masses, weights):
            fn = 'ana_datamc_dy%i%s.root' % (m, fset) if m != 20 else 'ana_datamc_zmumu.root'
            fn = hists_dir + fn
            f = ROOT.TFile(fn)
            d = f.Our2012MuonsPlusMuonsMinusHistos

            h = d.Get('DileptonMass').Clone('dy%i%s' % (m, fset))
            h.Rebin(rebin)
            h.GetXaxis().SetRangeUser(xmin, xmax)
            h.Scale(w)
            hists.append(h)

        htot = hists[0].Clone('htot')
        for j in xrange(1, len(hists)):
            htot.Add(hists[j])
        htot.SetTitle('')
        htots.append(htot)

    for i in xrange(0, len(htots)):
        if i == 0:
            htots[i].Draw()
        else:
            htots[i].Draw("same")
    htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
    htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    ps.save('mass_overlay')

    l = ROOT.TLine(xmin, 1., xmax, 1.)
    for i in xrange(1, len(htots)):
        ps.c.Clear()
        ratio = htots[i].Clone()
        ratio.Divide(htots[0])
        ratio.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
        ratio.GetYaxis().SetTitle('C%i / STARTUP' % i)
        ratio.SetMinimum(0.5)
        ratio.SetMaximum(1.5)
        ratio.SetMarkerStyle(20)
        ratio.Draw("pe")
        ffunc = ROOT.TF1('ffunc', '[0] + x*[1] + x*x*[2]', xmin, xmax)
        ratio.Fit(ffunc, 'IVR')
        l.Draw()
        ps.save('ratio%i_0' % i, log=False)
        if i > 1:
            ratio = htots[i].Clone()
            ratio.Divide(htots[i-1])
            ratio.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
            ratio.GetYaxis().SetTitle('C%i / C%i' % (i, i-1))
            ratio.SetMinimum(0.5)
            ratio.SetMaximum(1.5)
            ratio.SetMarkerStyle(24)
            ratio.SetMarkerColor(ROOT.kBlue)
            ffunc = ROOT.TF1('ffunc', '[0] + x*[1] + x*x*[2]', xmin, xmax)
            ffunc.SetLineColor(ROOT.kBlue)
            ratio.Fit(ffunc, 'IVR')
            ratio.Draw("pe")
            l.Draw()
            ps.save('ratio%i_%i' % (i, i-1), log=False)

#l = range(200, 2200, 200)
l = [200, 300, 400, 600, 800, 2400]
for lo in l:
    fit_it(lo, 2500)
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
#draw_overlay(fsets)
