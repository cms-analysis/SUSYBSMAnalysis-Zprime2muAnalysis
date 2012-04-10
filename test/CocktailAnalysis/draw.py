#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *

set_zp2mu_style()
rainbow_palette()

#ps = plot_saver('plots/cocktail')
f = ROOT.TFile('cocktail.root')

c = ROOT.TCanvas('c', '', 1100, 1400)
c.Print('cocktail.ps[')

for x in ['All', 'Barrel', 'Overlap', 'Endcap']:
    title = ROOT.TText(0.7, 0.7, x)
    
    d = f.Get('CocktailAnalyzer' + x)

    def do(s):
        h = d.Get(s)
        h.SetStats(0)
        h.Draw('hist text00')

    c.Clear()
    c.cd(0)
    title.Draw()
    c.Divide(1,3)
    c.cd(1)
    do('TMRCocktailChoice')
    c.cd(2)
    do('TunePCocktailChoice')
    c.cd(3)
    do('SigmaSwitchCocktailChoice')
    c.Print('cocktail.ps')

    def do(s):
        h = d.Get(s)
        h.GetXaxis().SetRangeUser(0, 4)
        h.GetYaxis().SetRangeUser(0, 4)
        h.Draw('colz')
        c.Update()
        st = h.FindObject('stats')
        st.SetX1NDC(0.63)
        st.SetX2NDC(0.88)
        st.Draw()
        return h

    c.Clear()
    c.cd(0)
    title.Draw()
    c.Divide(2,3)
    for i,x in enumerate(['GlobalVsTkOnly', 'GlobalVsTPFMS', 'GlobalVsPicky', 'TkOnlyVsTPFMS', 'TkOnlyVsPicky', 'TPFMSVsPicky']):
        c.cd(i+1)
        do('TrackLnChi2TailProb%s' % x)
    c.Print('cocktail.ps')
    
c.Print('cocktail.ps]')

