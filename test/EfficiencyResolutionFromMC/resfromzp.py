#!/usr/bin/env python

# (py resfromzp.py >! plots/out.resfromzp) && tlp plots/*resfromzp

import os, sys
from array import array
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)

ps = plot_saver('plots/resfromzp')

res_vs_mass = []
def drawMassReso(basein='Resolutioninner',refit='tkonly',space='DileptonMassRes',nrms=1.5):
    sigma = []
    esigma = []

    masses = [750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]
    chnam, val, err, xlolim, xhilim, iuint = ROOT.TString(''), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Long(1)
    baseout = "%s_%s"%(refit,space)
    for m in masses:
        f = ROOT.TFile('crab/effres_histos/Mu40/ana_effres_zp%i.root' % m)
        h = f.Get(basein).Get(space)
        fcn = ROOT.TF1('mg', 'gaus', h.GetMean() - nrms*h.GetRMS(), h.GetMean() + nrms*h.GetRMS() )
        h.Fit(fcn, 'ILVR')
        ps.save('%s_%i' % (baseout, m))
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

    print 'Fitting resolution-vs-mass curve:'
    g = ROOT.TGraphErrors(l, masses, sigma, emasses, esigma)
    # Define pol2 explicitly instead of using "pol2" to constrain one of
    # its parameters and prevent fit from exploding.
    poly2 = ROOT.TF1('poly2', '[0] + x*[1] + x*x*[2]')
    poly2.SetParLimits(0, 0., 1.);
    g.Fit(poly2)
    g.GetXaxis().SetTitle('dilepton mass (GeV)')
    g.GetYaxis().SetTitle('dilepton relative mass resolution')
    g.GetYaxis().SetTitleOffset(1.2);
    g.SetTitle('')
    g.Draw("AP")
    ps.c.Update()
    s = g.GetListOfFunctions().FindObject("stats")
    s.SetOptFit(111)
    s.SetY1NDC(0.15)
    s.SetY2NDC(0.35)
    s.SetX1NDC(0.50)
    s.SetX2NDC(0.90)
    ps.c.Update()
    ps.save(baseout, log=False)
    if space == 'DileptonMassRes':
        res_vs_mass.append((refit, g))

for basein, baseout in [('Resolutioninner', 'tkonly'), ('Resolutiontunep', 'tunep'), ('Resolutiontunepnew', 'tunep_new')]:
    for histname in ['DileptonMassRes','DileptonMassInvRes']:
        if 'Inv' in histname: drawMassReso(basein,baseout,histname,3.0)
        else: drawMassReso(basein,baseout,histname)

#- Compare the resolution given by different methods
ps.c.Clear()
# ROOT.gStyle.SetOptFit(0)
lg = ROOT.TLegend(0.50, 0.15, 0.90, 0.35)
for refit, graph in res_vs_mass:
#    s = graph.GetListOfFunctions().FindObject("stats")
#    s.SetOptFit(0)
    graph.SetMaximum(0.17)
    graph.SetMarkerSize(1.0)
    fcn = graph.GetFunction('poly2')
    if refit == 'tkonly':
        graph.SetMarkerStyle(20)
        graph.SetLineColor(ROOT.kRed)
        graph.SetMarkerColor(graph.GetLineColor())
        fcn.SetLineColor(graph.GetLineColor())
        graph.Draw("AP")
        lg.AddEntry(graph, 'Tracker-only', 'LPE')
    if refit == 'tunep':
        graph.SetMarkerStyle(24)
        graph.SetLineColor(ROOT.kBlue)
        graph.SetMarkerColor(graph.GetLineColor())
        fcn.SetLineColor(graph.GetLineColor())
        graph.Draw("P same")
        lg.AddEntry(graph, 'Tune P', 'LPE')
    if refit == 'tunep_new':
        graph.SetMarkerStyle(21)
        graph.SetLineColor(ROOT.kGreen)
        graph.SetMarkerColor(graph.GetLineColor())
        fcn.SetLineColor(graph.GetLineColor())
        graph.Draw("P same")
        lg.AddEntry(graph, 'New Tune P', 'LPE')
lg.Draw()
ps.save('res_vs_mass', log=False)
