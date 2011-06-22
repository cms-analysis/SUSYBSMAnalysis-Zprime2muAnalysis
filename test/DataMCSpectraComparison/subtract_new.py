#!/usr/bin/env python

import sys, os, glob
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetTitleX(0.12)
#ROOT.gStyle.SetTitleH(0.07)

ps = plot_saver('plots/subtract_new_from_old_sel', pdf_log=True)

f = ROOT.TFile('ana_datamc_Run2011AMuonsOnly/ana_datamc_data.root')

hold = f.OurOldMuonsPlusMuonsMinusHistos.Get('DileptonMass')
hnew = f.OurNewMuonsPlusMuonsMinusHistos.Get('DileptonMass')

for h in (hold,hnew):
    h.Rebin(5)

hdiff = hold.Clone('hdiff')
hdiff.Add(hnew, -1)

hdiff.SetTitle(';m(#mu^{+}#mu^{-}) [GeV];Events/5 GeV')
hdiff.GetXaxis().SetRangeUser(60, 800)
hdiff.SetMarkerStyle(20)
hdiff.SetMarkerSize(0.8)
hdiff.SetStats(0)
hdiff.Draw('e')
ps.save('MuonsPlusMuonsMinus_OldSelMinusNew')
