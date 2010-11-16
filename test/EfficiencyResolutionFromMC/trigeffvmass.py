#!/usr/bin/env python

import sys, os

samples = ['dy20', 'dy120', 'dy200', 'dy500', 'dy800', 'zp500', 'zp750', 'zp1000', 'zp1250', 'zp1500', 'zp1750']
dy = [(x, int(x.replace('dy',''))) for x in samples if 'dy' in x]
zp = [(x, int(x.replace('zp',''))) for x in samples if 'zp' in x]
if 'vbtf' in sys.argv:
    which = 'VBTFEfficiencyFromMC'
    plot_dir = 'plots/trigeffvsmassmctruth_vbtf'
else:
    which = 'EfficiencyFromMC'
    plot_dir = 'plots/trigeffvsmassmctruth'

from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)

types = ['Acceptance', 'RecoWrtAcc', 'RecoWrtAccTrig', 'TotalReco', 'L1OrEff', 'HLTOrEff', 'TotalTrigEff']
rebin_factor = 100
make_individual_effs = False

files = {}
plots = {}

#samples, types = samples[0:1], types[0:1]

samples_totals = []

print '%30s%30s%30s%30s%30s' % ('sample', 'type', 'num', 'den', 'eff')
for sample in samples:
    f = files[sample] = ROOT.TFile('ana_effres_%s.root' % sample)
    d = f.Get(which)
    ps = plot_saver(os.path.join(plot_dir, sample))

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
        samples_totals.append((sample, t, cnum, cden))
        sys.stdout.flush()

ps = plot_saver(plot_dir)

totals_histos = {}
for sample, t, cnum, cden in samples_totals:
    if 'zp' not in sample and sample != 'dy20':
        continue
    if not totals_histos.has_key(t):
        totals_histos[t] = ROOT.TH1F(t + '_totals_num', '', 2000, 0, 2000), ROOT.TH1F(t + '_totals_den', '', 2000, 0, 2000)
    hnum, hden = totals_histos[t]
    if 'zp' in sample:
        mass = int(sample.replace('zp', ''))
    else:
        mass = 90
    mass = hden.FindBin(mass)
    hnum.SetBinContent(mass, cnum)
    hden.SetBinContent(mass, cden)

for t in sorted(totals_histos.iterkeys()):
    eff = binomial_divide(*totals_histos[t])
    if 'Accept' in t:
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
L1Or.GetXaxis().SetRangeUser(100, 2000)
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
lg = ROOT.TLegend(0.42, 0.13, 0.92, 0.33)
lg.AddEntry(RecoWrtAccTrig, 'wrt triggered events in acceptance', 'LPE')
lg.AddEntry(RecoWrtAcc, 'wrt events in acceptance', 'LPE')
lg.AddEntry(TotalReco, 'total', 'LPE')
RecoWrtAccTrig.GetXaxis().SetRangeUser(100, 2000)
RecoWrtAccTrig.GetYaxis().SetRangeUser(0.5, 1.01)
RecoWrtAccTrig.Draw('AP')
RecoWrtAcc.Draw('P same')
TotalReco.Draw('P same')
lg.Draw()
ps.save('summary_recoeff', log=False)
