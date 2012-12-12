#!/usr/bin/env python

# (py trigeffvmass.py >! plots/out.trigeffvsmassmctruth) && (py trigeffvmass.py vbtf >! plots/out.trigeffvsmassmctruth_vbtf) && tlp plots/*trigeffvsmassmctruth*

import sys, os
from array import array

samples = ['dy60', 'dy120', 'dy200', 'dy500', 'dy800', 'dy1000', 'dy1500', 'dy2000', 'zp750', 'zp1000', 'zp1250', 'zp1500', 'zp1750', 'zp2000', 'zp2250', 'zp2500', 'zp2750', 'zp3000']

kind = [x for x in sys.argv[1:] if os.path.isdir(x)]

dy = [(x, int(x.replace('dy','').replace('_c1',''))) for x in samples if 'dy' in x]
zp = [(x, int(x.replace('zp','').replace('_c1',''))) for x in samples if 'zp' in x]
if 'vbtf' in sys.argv:
    which = 'VBTFEfficiencyFromMC'
    plot_dir = 'plots/trigeffvsmassmctruth_vbtf'
else:
    which = 'EfficiencyFromMC'
    plot_dir = 'plots/trigeffvsmassmctruth'

if kind:
    kind = kind[0]
    plot_dir += '_' + kind
else:
    kind = 'normal'

from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)

types = ['AccNoPt', 'Acceptance', 'RecoWrtAcc', 'RecoWrtAccTrig', 'TotalReco', 'L1OrEff', 'HLTOrEff', 'TotalTrigEff']
rebin_factor = 50
make_individual_plots = False
make_individual_effs = False

files = {}
plots = {}

#samples, types = samples[0:1], types[0:1]

samples_totals = []

print '%30s%30s%30s%30s%30s' % ('sample', 'type', 'num', 'den', 'eff')
for sample in samples:
    f = files[sample] = ROOT.TFile(os.path.join(kind, 'ana_effres_%s.root' % sample))
    d = f.Get(which)
    if make_individual_plots or make_individual_effs:
        ps = plot_saver(os.path.join(plot_dir, sample), pdf=True)

    for t in types:
        num = d.Get('Num' + t)
        den = d.Get('Den' + t)

        num.Rebin(rebin_factor)
        den.Rebin(rebin_factor)

        if make_individual_plots:
            num.Draw()
            ps.save(t + 'Num')
            den.Draw()
            ps.save(t + 'Den')

            num.SetLineColor(ROOT.kRed)
            den.SetLineColor(ROOT.kBlack)
            den.Draw()
            num.Draw('same')
            ps.save(t + 'Both')

        if make_individual_effs:
            eff = binomial_divide(num, den)
            eff.Draw('AP')
            if t not in ['AccNoPt', 'Acceptance']:
                eff.GetYaxis().SetRangeUser(0.85, 1.01)
            else:
                eff.GetYaxis().SetRangeUser(0.5, 1.05)
            ps.save(t, log=False)
        else:
            eff = None

        plots[(sample, t)] = (num, den, eff)

        nb = num.GetNbinsX()
        cnum = num.Integral(0, nb+1)
        cden = den.Integral(0, nb+1)

        print '%30s%30s%30f%30f%30f' % (sample, t, cnum, cden, cnum/cden)
        samples_totals.append((sample, t, cnum, cden))
        sys.stdout.flush()

ps = plot_saver(plot_dir, pdf=True)

totals_histos = {}
for sample, t, cnum, cden in samples_totals:
    if 'zp' not in sample and sample != 'dy60':
        continue
    if not totals_histos.has_key(t):
        totals_histos[t] = ROOT.TH1F(t + '_totals_num', '', 3100, 0, 3100), ROOT.TH1F(t + '_totals_den', '', 3100, 0, 3100)
    hnum, hden = totals_histos[t]
    if 'zp' in sample:
        mass = int(sample.replace('zp', '').replace('_c1',''))
    else:
        mass = 90
    mass = hden.FindBin(mass)
    hnum.SetBinContent(mass, cnum)
    hden.SetBinContent(mass, cden)

for t in sorted(totals_histos.iterkeys()):
    hnum, hden = totals_histos[t]
    eff = binomial_divide(hnum, hden)
    totals_histos[t + '_eff'] = eff
    if t.startswith('Acc'):
        eff.GetYaxis().SetRangeUser(0.4, 1.01)
    elif 'Reco' in t:
        eff.GetYaxis().SetRangeUser(0., 1.01)
    else:
        eff.GetYaxis().SetRangeUser(0.75, 1.01)
    eff.Draw('AP')
    ps.save(t + '_totals', log=False)

def do_overlay(name, samples, which):
    if not samples:
        return
    x = 'AP'
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, 8, 40]
    for i, ((sample, mass), color) in enumerate(zip(samples, colors)):
        g = plots[(sample, which)][-1]
        g.SetLineColor(color)
        g.Draw(x)
        g.GetXaxis().SetLimits(40, 2500)
        if which not in ['AccNoPt', 'Acceptance']:
            g.GetYaxis().SetRangeUser(0.85, 1.05)
        else:
            g.GetYaxis().SetRangeUser(0.4, 1.05)
        x = 'P same'
    ps.save('%s_%s_overlay' % (name, which))

if make_individual_effs:
    do_overlay('dy',dy,'Acceptance')
    do_overlay('zp',zp,'Acceptance')
    do_overlay('dy',dy,'TotalTrigEff')
    do_overlay('zp',zp,'TotalTrigEff')

def do_replace_or_add(name, samples, which, add=False):
    if not samples:
        return
    num = plots[(samples[0][0], which)][0].Clone('%s_numtot' % name)
    den = plots[(samples[0][0], which)][1].Clone('%s_dentot' % name)
    assert(num.GetNbinsX() == den.GetNbinsX())
    for i in xrange(1, len(samples)):
        n = plots[(samples[i][0], which)][0]
        d = plots[(samples[i][0], which)][1]
        assert(n.GetNbinsX() == d.GetNbinsX() and n.GetNbinsX() == num.GetNbinsX())
        for j in xrange(0, n.GetNbinsX()+1):
            if add:
                num.SetBinContent(j, num.GetBinContent(j) + n.GetBinContent(j))
                den.SetBinContent(j, den.GetBinContent(j) + d.GetBinContent(j))
            else: # replace
                if d.GetBinContent(j):
                    num.SetBinContent(j, n.GetBinContent(j))
                    den.SetBinContent(j, d.GetBinContent(j))

    num.SetLineColor(ROOT.kRed)
    den.SetLineColor(ROOT.kBlack)
    den.Draw('hist')
    num.Draw('hist same')
    ps.save('%s_%s_%s_both' % (name, which, 'add' if add else 'replace'))

    eff = binomial_divide(num, den)
    eff.Draw('AP')
    eff.GetXaxis().SetLimits(40, 2500)
    if which.startswith('Acc'):
        eff.GetYaxis().SetRangeUser(0.2, 1.01)
    elif 'Reco' in which:
        eff.GetYaxis().SetRangeUser(0., 1.01)
    else:
        eff.GetYaxis().SetRangeUser(0.75, 1.01)
        
    eff.SetTitle('')
    eff.GetXaxis().SetTitle('dimuon invariant mass (GeV)')
    eff.GetYaxis().SetLabelSize(0.03)
    eff.GetYaxis().SetTitle('acceptance' if which in ['AccNoPt', 'Acceptance'] else 'efficiency')
    ps.save('%s_%s_%s' % (name, which, 'add' if add else 'replace'), log=False)
    return eff

final_plots = {}
for t in types:
    final_plots[t] = do_replace_or_add('dy', dy, t, add=True)

for x,y in final_plots.iteritems():
    exec '%s = y' % x
   
L1Or = final_plots['L1OrEff']
HLTOr = final_plots['HLTOrEff']
Total = final_plots['TotalTrigEff']

'''
L1Single = final_plots['L1Path_0_L1_SingleMu7']
L1Single.SetLineColor(ROOT.kRed)
L1Or.SetLineColor(ROOT.kBlue)
lg = ROOT.TLegend(0.13, 0.13, 0.44, 0.30)
lg.AddEntry(L1Single, 'L1 single', 'LE')
lg.AddEntry(L1Or, 'L1 single OR double', 'LE')
L1Single.Draw('AP')
L1Or.Draw('P same')
lg.Draw()
ps.save('summary_l1', log=False)

HLTSingle = final_plots['HLTPath_0_HLT_Mu11']
HLTSingle.SetLineColor(ROOT.kRed)
HLTOr.SetLineColor(ROOT.kBlue)
HLTSingle.Draw('AP')
HLTOr.Draw('P same')
lg = ROOT.TLegend(0.13, 0.13, 0.44, 0.30)
lg.AddEntry(HLTSingle, 'HLT single', 'LE')
lg.AddEntry(HLTOr, 'HLT single OR double', 'LE')
lg.Draw()
ps.save('summary_hlt', log=False)
'''

for h,c,m in [(L1Or, ROOT.kRed, 20), (HLTOr, ROOT.kGreen+2, 21), (Total, ROOT.kBlue, 22)]:
    h.SetMarkerStyle(m)
    h.SetMarkerColor(c)
    h.SetMarkerSize(1.2)
    h.SetLineColor(c)
lg = ROOT.TLegend(0.13, 0.13, 0.44, 0.30)
lg.AddEntry(L1Or, 'L1', 'LPE')
lg.AddEntry(HLTOr, 'HLT', 'LPE')
lg.AddEntry(Total, 'L1+HLT', 'LPE')
L1Or.GetXaxis().SetRangeUser(100, 2500)
L1Or.GetYaxis().SetRangeUser(0.93, 1.005)
L1Or.Draw('AP')
HLTOr.Draw('P same')
Total.Draw('P same')
lg.Draw()
ps.save('summary_total', log=False)

for h,c,m in [(RecoWrtAccTrig, ROOT.kRed, 20), (RecoWrtAcc, ROOT.kGreen+2, 21), (TotalReco, ROOT.kBlue, 22)]:
    h.SetMarkerStyle(m)
    h.SetMarkerColor(c)
    h.SetMarkerSize(1.2)
    h.SetLineColor(c)
lg = ROOT.TLegend(0.20, 0.13, 0.92, 0.33)
lg.AddEntry(RecoWrtAccTrig, 'wrt triggered events in acceptance', 'LPE')
lg.AddEntry(RecoWrtAcc, 'wrt events in acceptance', 'LPE')
lg.AddEntry(TotalReco, 'total', 'LPE')
RecoWrtAccTrig.GetXaxis().SetRangeUser(0, 2500)
RecoWrtAccTrig.GetYaxis().SetRangeUser(0.1, 1.01)
RecoWrtAccTrig.Draw('AP')
RecoWrtAcc.Draw('P same')
TotalReco.Draw('P same')
fitwindow = 200,2500
fcn = ROOT.TF1('fcn', '[0] + [1]/(x + [2])**3', *fitwindow)
fcn.SetParNames("a", "b", "c")
fcn.SetLineColor(TotalReco.GetLineColor())
lg.AddEntry(fcn, 'fit to a + b/pow(m+c, 3) for m in (200,2500)', 'L')
lg.Draw()
ROOT.gStyle.SetOptFit(1111)
TotalReco.Fit(fcn, 'VR')
ps.c.Update()
s = TotalReco.GetListOfFunctions().FindObject('stats')
s.SetFitFormat('5.3g')
s.SetX1NDC(0.49)
s.SetY1NDC(0.35)
s.SetX2NDC(0.92)
s.SetY2NDC(0.59)
ps.save('summary_recoeff', log=False)

residuals = TotalReco.Clone('residuals')
x,y = ROOT.Double(), ROOT.Double()
for i in xrange(TotalReco.GetN()):
    TotalReco.GetPoint(i, x, y)
    f = fcn.Eval(x)
    residuals.SetPoint(i, x, f/y-1)
residuals.GetXaxis().SetRangeUser(0, 2500)
residuals.SetTitle(';dimuon invariant mass (GeV);relative residual f/h-1')
residuals.Draw('AP')
residuals.Fit('pol1', 'VR', '', *fitwindow)
fcn = residuals.GetFunction('pol1')
fcn.SetParNames('a', 'b')
ps.c.Update()
residuals.GetYaxis().SetRangeUser(-0.1, 0.1)
residuals.GetYaxis().SetTitleOffset(1.2)
lg = ROOT.TLegend(0.14, 0.32, 0.44, 0.39)
lg.AddEntry(fcn, 'Fit to a + bx', 'L')
lg.Draw()
s = residuals.GetListOfFunctions().FindObject('stats')
s.SetX1NDC(0.14)
s.SetY1NDC(0.14)
s.SetX2NDC(0.54)
s.SetY2NDC(0.31)
ps.save('totalreco_fit_residuals', log=False)

# Dump the values of the total reconstruction curve.  Also print effs.
# for 60-120 and 120-200 GeV bins from the total counts.
print '\nTotal efficiencies (acceptance times trigger+reconstruction+selection efficiencies)'
print '%20s%10s%20s' % ('mass range', 'eff', '68%CL interval')
x,y = ROOT.Double(), ROOT.Double()
g = totals_histos['TotalReco_eff']
g.GetPoint(0, x, y)
print '%20s%10.4f%20s' % ( '60-120', y, '[%5.4f, %5.4f]' % (y-g.GetErrorYlow(0), y+g.GetErrorYhigh(0)))
g = TotalReco
n = g.GetN()
for i in xrange(n):
    g.GetPoint(i, x, y)
    exl = g.GetErrorXlow(i)
    exh = g.GetErrorXhigh(i)
    eyl = g.GetErrorYlow(i)
    eyh = g.GetErrorYhigh(i)
    print '%20s%10.4f%20s' % ('%.f-%.f' % (x-exl, x+exh), y, '[%5.4f, %5.4f]' % (y-eyl, y+eyh))
