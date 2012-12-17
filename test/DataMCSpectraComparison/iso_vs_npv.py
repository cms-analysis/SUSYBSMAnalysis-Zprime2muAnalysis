#!/usr/bin/env python

from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ROOT.gStyle.SetOptStat(110010)

ps = plot_saver('plots/iso_vs_npv', pdf=True)

f = ROOT.TFile('data/Run2012MuonsOnly/ana_datamc_data.root')
t = f.SimpleNtupler.Get('t')

max_npv = 30
nbins = 6
per = max_npv/nbins

hs = []

t.SetAlias('OurSelNoIso', 'loose_no_iso_0 && loose_no_iso_1 && triggerMatched && extraDimuonCuts && GoodData && OppSign')
t.Draw('nvertices>>hnvertices(30, 0, 30)', 'OurSelNoIso && dil_mass > 60 && dil_mass < 120')
ps.save('nvertices')
t.Draw('nvertices>>hnvertices2(6, 0, 30)', 'OurSelNoIso && dil_mass > 60 && dil_mass < 120')
ROOT.hnvertices2.Draw('histo text00')
ps.save('nvertices_rebin')

for tkonly in [0, 1]:
    print "apples"
    if tkonly:
        t.SetAlias('iso', 'lep_sumPt/lep_tk_pt')
        cut = 0.1
        which = 'tkonly'
    else:
        t.SetAlias('iso', '(lep_sumPt+lep_emEt+lep_hadEt)/lep_tk_pt')
        cut = 0.15
        which = 'tkpluscalo'

    t.Draw('iso:nvertices>>hscatter%i(%i,0,%i,100,0,5)' % (tkonly, nbins, max_npv), 'OurSelNoIso && dil_mass > 60 && dil_mass < 120')
    h = getattr(ROOT, 'hscatter%i' % tkonly)
    h.SetTitle('selected muons from Z#rightarrow#mu#mu (60 < m(#mu#mu) < 120 GeV);number of reconstructed primary vertices;%s relative isolation' % ('tracker-only' if tkonly else 'tracker-plus-calo.'))
    for i in xrange(nbins):
        h.GetXaxis().SetBinLabel(i+1, '%i-%i' % (i*per, (i+1)*per))
    h.Draw('box')
    ps.save('scatter_%s' % which)

    habove = ROOT.TH1F('habove%i' % tkonly, '', nbins, 0, max_npv)
    htotal = ROOT.TH1F('htotal%i' % tkonly, '', nbins, 0, max_npv)
    for i in xrange(nbins):
        proj = h.ProjectionY('hproj%i' % i, i+1, i+1)
        habove.SetBinContent(i+1, proj.Integral(proj.FindBin(cut), proj.GetNbinsX()+2))
        htotal.SetBinContent(i+1, proj.GetEntries())

    hfrac = binomial_divide(habove, htotal)
    hfrac.SetTitle('selected muons from Z#rightarrow#mu#mu (60 < m(#mu#mu) < 120 GeV);number of reconstructed primary vertices;fraction of events with rel. iso. > %.2f' % cut)
    bnr = hfrac.GetXaxis().GetNbins()/(hfrac.GetN()+0.5)
    for i in xrange(nbins):
        hfrac.GetXaxis().SetBinLabel(int((i+0.55)*bnr), '%i-%i' % (i*per, (i+1)*per))
    hfrac.GetXaxis().LabelsOption('h')
    hfrac.GetYaxis().SetTitleOffset(1.2)
    hfrac.GetYaxis().SetLabelSize(0.025)
    hfrac.Draw('ap')
    ps.save('%s_iso_above_%s_vs_npv_frac' % (which, str(cut).replace('.','')))

    hs.append(hfrac.Clone(which))

tkpcal, tkonly = hs
tkpcal.SetLineColor(ROOT.kRed)
tkpcal.SetMarkerColor(ROOT.kRed)
tkpcal.SetLineWidth(2)
tkonly.SetLineWidth(2)
tkpcal.SetMarkerStyle(21)
tkonly.SetMarkerStyle(22)
tkpcal.SetMarkerSize(1.1)
tkonly.SetMarkerSize(1.2)
tkpcal.Draw('APEZ')
tkpcal.GetYaxis().SetTitle('fraction of events failing the rel. iso. cut')
tkpcal.GetYaxis().SetRangeUser(0.005, 0.05)
tkonly.Draw('PEZ')
l = ROOT.TLegend(0.15, 0.74, 0.70, 0.86)
l.AddEntry(tkpcal, 'tracker-plus-calorimeters, cut = 0.15', 'LPE')
l.AddEntry(tkonly, 'tracker-only, cut = 0.1', 'LPE')
l.Draw()
ps.save('overlay')
