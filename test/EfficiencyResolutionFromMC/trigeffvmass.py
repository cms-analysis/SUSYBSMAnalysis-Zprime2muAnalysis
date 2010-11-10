#!/usr/bin/env python

import sys, os

samples = ['zmumu', 'dy120', 'dy200', 'dy500', 'dy800', 'zp1000', 'zp1250', 'zp1500', 'zp1750']
dy = [('zmumu', 40)] + [(x, int(x.replace('dy',''))) for x in samples if 'dy' in x] # and 'dy120' not in x]
zp = [(x, int(x.replace('zp',''))) for x in samples if 'zp' in x]

from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)

#types = ['Acceptance', 'L1OrEff', 'L1AndEff', 'HLTOrEff', 'HLTAndEff', 'TotalTrigEff']
types = ['Acceptance', 'RecoWrtAcc', 'RecoWrtAccTrig', 'TotalReco', 'L1Path_0_L1_SingleMu7', 'L1Path_1_L1_DoubleMu3', 'L1OrEff', 'HLTPath_0_HLT_Mu9', 'HLTPath_1_HLT_DoubleMu3', 'HLTOrEff', 'TotalTrigEff']
rebin_factor = 100
make_individual_effs = False

files = {}
plots = {}

#samples, types = samples[0:1], types[0:1]

print '%30s%30s%30s%30s%30s' % ('sample', 'type', 'num', 'den', 'eff')
for sample in samples:
    f = files[sample] = ROOT.TFile('ana_effres_%s.root' % sample)
    d = f.EfficiencyFromMC
    ps = plot_saver('plots/trigeffvsmassmctruth/%s' % sample)

    for t in types:
        num = d.Get('Num' + t)
        den = d.Get('Den' + t)

        num.Rebin(rebin_factor)
        den.Rebin(rebin_factor)
                
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
            if t != 'Acceptance':
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
        sys.stdout.flush()

ps = plot_saver('plots/trigeffvsmassmctruth')

def do_overlay(name, samples, which):
    x = 'AP'
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, 8, 40]
    for i, ((sample, mass), color) in enumerate(zip(samples, colors)):
        g = plots[(sample, which)][-1]
        g.SetLineColor(color)
        g.Draw(x)
        g.GetXaxis().SetLimits(40, 2000)
        if which != 'Acceptance':
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
    den.Draw()
    num.Draw('same')
    ps.save('%s_%s_%s_both' % (name, which, 'add' if add else 'replace'))

    eff = binomial_divide(num, den)
    eff.Draw('AP')
    eff.GetXaxis().SetLimits(40, 2000)
    if 'Accept' in which:
        eff.GetYaxis().SetRangeUser(0.4, 1.01)
    elif 'Reco' in which:
        eff.GetYaxis().SetRangeUser(0., 1.01)
    else:
        eff.GetYaxis().SetRangeUser(0.75, 1.01)
        
    eff.SetTitle('')
    eff.GetXaxis().SetTitle('dimuon invariant mass (GeV)')
    eff.GetYaxis().SetLabelSize(0.03)
    eff.GetYaxis().SetTitle('acceptance' if which == 'Acceptance' else 'efficiency')
    ps.save('%s_%s_%s' % (name, which, 'add' if add else 'replace'), log=False)
    return eff

if False:
    do_replace_or_add('dy', dy, 'Acceptance')
    do_replace_or_add('dy', dy, 'TotalTrigEff')
    do_replace_or_add('zp', zp, 'Acceptance')
    do_replace_or_add('zp', zp, 'TotalTrigEff')

final_plots = {}
for t in types:
    final_plots[t] = do_replace_or_add('dy', dy, t, add=True)
    do_replace_or_add('zp', zp, t, add=True)

for x,y in final_plots.iteritems():
    exec '%s = y' % x
   
L1Single = final_plots['L1Path_0_L1_SingleMu7']
L1Double = final_plots['L1Path_1_L1_DoubleMu3']
L1Or = final_plots['L1OrEff']
HLTSingle = final_plots['HLTPath_0_HLT_Mu9']
HLTDouble = final_plots['HLTPath_1_HLT_DoubleMu3']
HLTOr = final_plots['HLTOrEff']
Total = final_plots['TotalTrigEff']
   
L1Single.SetLineColor(ROOT.kRed)
L1Double.SetLineColor(ROOT.kGreen+2)
L1Or.SetLineColor(ROOT.kBlue)
lg = ROOT.TLegend(0.13, 0.13, 0.44, 0.30)
lg.AddEntry(L1Single, 'L1 single', 'LE')
lg.AddEntry(L1Double, 'L1 double', 'LE')
lg.AddEntry(L1Or, 'L1 single OR double', 'LE')
L1Single.Draw('AP')
L1Double.Draw('P same')
L1Or.Draw('P same')
lg.Draw()
ps.save('summary_l1', log=False)

HLTSingle.SetLineColor(ROOT.kRed)
HLTDouble.SetLineColor(ROOT.kGreen+2)
HLTOr.SetLineColor(ROOT.kBlue)
HLTSingle.Draw('AP')
HLTDouble.Draw('P same')
HLTOr.Draw('P same')
lg = ROOT.TLegend(0.13, 0.13, 0.44, 0.30)
lg.AddEntry(HLTSingle, 'HLT single', 'LE')
lg.AddEntry(HLTDouble, 'HLT double', 'LE')
lg.AddEntry(HLTOr, 'HLT single OR double', 'LE')
lg.Draw()
ps.save('summary_hlt', log=False)

L1Or.SetLineColor(ROOT.kRed)
HLTOr.SetLineColor(ROOT.kGreen+2)
Total.SetLineColor(ROOT.kBlue)
lg = ROOT.TLegend(0.13, 0.13, 0.44, 0.30)
lg.AddEntry(L1Or, 'L1', 'LE')
lg.AddEntry(HLTOr, 'HLT', 'LE')
lg.AddEntry(Total, 'L1+HLT', 'LE')
L1Or.Draw('AP')
HLTOr.Draw('P same')
Total.Draw('P same')
lg.Draw()
ps.save('summary_total', log=False)

RecoWrtAccTrig.SetLineColor(ROOT.kRed)
RecoWrtAcc.SetLineColor(ROOT.kGreen+2)
Total.SetLineColor(ROOT.kBlue)
lg = ROOT.TLegend(0.42, 0.13, 0.92, 0.33)
lg.AddEntry(RecoWrtAccTrig, 'wrt triggered events in acceptance', 'LE')
lg.AddEntry(RecoWrtAcc, 'wrt events in acceptance', 'LE')
lg.AddEntry(TotalReco, 'total', 'LE')
RecoWrtAccTrig.Draw('AP')
RecoWrtAcc.Draw('P same')
TotalReco.Draw('P same')
lg.Draw()
ps.save('summary_recoeff', log=False)
