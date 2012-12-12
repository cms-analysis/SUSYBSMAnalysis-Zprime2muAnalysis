#!/usr/bin/env python

# (py resfromdy.py >! plots/out.resfromdy) && tlp plots/*resfromdy

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)

ps = plot_saver('plots/resfromdy')

f = ROOT.TFile('dyall_c1.root')
d = f.Resolutiontunepnew
h = d.Get('DileptonMassResVMass')
h = make_rms_hist(h)

h.SetTitle('')
h.GetXaxis().SetTitle('dilepton mass (GeV)')
h.GetYaxis().SetTitle('dilepton relative mass resolution')
h.GetYaxis().SetLabelSize(0.02)
h.GetXaxis().SetRangeUser(120, 2500)
h.Draw()

line = ROOT.TF1('line', '[0] + x*[1] + x*x*[2]', 120, 2500)
h.Fit(line, 'IVR')
ps.c.Update()
s = h.GetListOfFunctions().FindObject("stats")
s.SetOptStat(0)
s.SetOptFit(111)
s.SetY1NDC(0.15)
s.SetY2NDC(0.40)
s.SetX1NDC(0.50)
s.SetX2NDC(0.95)
ps.c.Update()

if 1:
    chnam, val, err, xlolim, xhilim, iuint = ROOT.TString(''), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Long(1)

    z = [0,0,0]
    ROOT.gMinuit.mnpout(0, chnam, val, err, xlolim, xhilim, iuint)
    z[0] = float(val)
    ROOT.gMinuit.mnpout(1, chnam, val, err, xlolim, xhilim, iuint)
    z[1] = float(val)
    ROOT.gMinuit.mnpout(2, chnam, val, err, xlolim, xhilim, iuint)
    z[2] = float(val)
    print z

    t = ROOT.TLatex(0.15, 0.8, '#sigma(M)/M = %.2e + %.2e M + %.2e M  ^{2}' % tuple(z))
    t.SetNDC()
    t.Draw()

ps.save('mass_res')
