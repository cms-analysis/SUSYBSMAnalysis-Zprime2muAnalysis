#!/usr/bin/env python

import sys, os, glob
from collections import defaultdict
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import cumulative_histogram, get_integral, move_above_into_bin, plot_saver, real_hist_max, real_hist_min, set_zp2mu_style, ROOT
set_zp2mu_style()
ROOT.gStyle.SetLineWidth(2)

from samples import *

rebin_factor = 5
x_axis_limits = 40, 1000
x_axis_limits2 = 40, 500
to_compare = 'DileptonMass'
global_rescale = 3273/3404.6 if False else None
draw_zssm = True

joins = [(s.name, 'jets') for s in samples if 'qcd' in s.name]
joins += [(x, 'jets') for x in ['inclmu15', 'wmunu', 'wjets']]
joins += [(x, 't#bar{t}-like') for x in ['singletop_tW', 'ztautau', 'ww', 'wz', 'zz']]
joins += [(s.name, 'Z #rightarrow #mu^{+}#mu^{-}') for s in samples if 'dy' in s.name]
joins += [('zmumu', 'Z #rightarrow #mu^{+}#mu^{-}')]

joins = dict(joins)
joins_colors = {'jets': 801, 't#bar{t}-like': 4, 'Z #rightarrow #mu^{+}#mu^{-}': 7}

histo_dir = None
for x in sys.argv:
    if 'xax' in x:
        xax0, xax1 = x.split('=')[1].split(',')
        x_axis_limits = int(xax0), int(xax1)
    elif 'compare' in x:
        to_compare = x.split('=')[1]
    elif os.path.isdir(x):
        histo_dir = x
assert(histo_dir is not None)

data_fns = glob.glob(os.path.join(histo_dir, 'ana_datamc_data*.root'))
if len(data_fns) == 1:
    data_fn = data_fns[0]
    hadd_tmp = False
else:
    # just hadd to tmp file for now, easier
    data_fn = 'anadatamcdatahaddtmp.root'
    hadd_tmp = True
    os.system('hadd -f %s %s' % (data_fn, ' '.join(data_fns)))
fdata = ROOT.TFile(data_fn)

def parse_lumi_from_log(fn):
    # Yay for fragile parsing!
    this = False
    for line in open(fn):
        if this:
            x = float(line.split()[-2])/1e6
            print fn, x
            return x
        if line == '-------------------------------------------------------------------\n':
            this = True

int_lumi = sum(parse_lumi_from_log(x.replace('.root', '.lumi')) for x in data_fns)
print 'total lumi from data: %.1f/pb' % int_lumi

subtitleize= {
    'MuonsPlusMuonsMinus': '#mu^{+}#mu^{-}',
    'MuonsPlusMuonsPlus':  '#mu^{+}#mu^{+}',
    'MuonsMinusMuonsMinus': '#mu^{-}#mu^{-}',
    'MuonsSameSign': '#mu^{#pm}#mu^{#pm}',
    'ElectronsPlusElectronsMinus': 'e^{+}e^{-}',
    'ElectronsPlusElectronsPlus': 'e^{+}e^{+}',
    'ElectronsMinusElectronsMinus': 'e^{-}e^{-}',
    'ElectronsSameSign': 'e^{#pm}e^{#pm}',
    'MuonsPlusElectronsMinus': '#mu^{+}e^{-}',
    'MuonsMinusElectronsPlus': '#mu^{-}e^{+}',
    'MuonsPlusElectronsPlus': '#mu^{+}e^{+}',
    'MuonsMinusElectronsMinus': '#mu^{-}e^{-}',
    'MuonsElectronsOppSign': '#mu^{#pm}e^{-+}',
    'MuonsElectronsSameSign': '#mu^{#pm}e^{#pm}',
    }
titleize = {
    'DileptonMass': 'm(%s) (%s)',
    'DileptonPt': '%s p_{T} (%s)'
    }
unitize = {
    'DileptonMass': 'GeV',
    'DileptonPt': 'GeV'
    }
yaxis = {
    'MuonsPlusMuonsMinus': (1e-3, None),
#    'MuonsSameSign': (5e-5, 2.5),
#    'MuonsElectronsOppSign': (2e-3, 40),
    }
use_yaxis = True

dileptons = ['MuonsPlusMuonsMinus', 'MuonsSameSign', 'MuonsElectronsOppSign', 'MuonsElectronsSameSign']
cutss = ['VBTF', 'Our', 'OurNoIso']

ROOT.TH1.AddDirectory(False)

def dir_name(c, d):
    return c + d + 'Histos'

sumMCfile = None # ROOT.TFile('sumMC.root', 'recreate')

pdir = 'plots/datamc'
if histo_dir != 'ana_datamc':
    pdir += '_' + histo_dir.replace('ana_datamc_', '')
if global_rescale is not None:
    pdir += '_normToZ'
ps = plot_saver(pdir)

if global_rescale is not None:
    for s in samples:
        s.partial_weight *= global_rescale

#samples = [s for s in samples if not s.name in ['ww', 'zz', 'wz', 'qcd500']]
for s in samples:
    s.histos = {}
    
for cuts in cutss:
    plot_dir = pdir + '/%s/%s' % (to_compare, cuts)
    ps.set_plot_dir(plot_dir)

    data = dict((d, getattr(fdata, dir_name(cuts, d)).Get(to_compare).Clone()) for d in dileptons)

    for dilepton in dileptons:
        if int_lumi > 39 and 'Electron' in dilepton:
            continue

        do_joins = True # ('Plus' in dilepton and 'Minus' in dilepton) or 'OppSign' in dilepton
        
        xax = x_axis_limits if dilepton == 'MuonsPlusMuonsMinus' else x_axis_limits2
        
        for sample in samples:
            # It would be more efficient to have the sample loop be
            # the outer one thanks to the file opening/closing, but
            # the code is cleaner this way.
            f = ROOT.TFile(os.path.join(histo_dir, 'ana_datamc_%s.root' % sample.name))
            sample.mass = getattr(f, dir_name(cuts, dilepton)).Get(to_compare).Clone()
            sample.mass.Rebin(rebin_factor)
            sample.mass.Scale(sample.partial_weight * int_lumi)
            move_above_into_bin(sample.mass, xax[1])
            sample.histos[cuts + dilepton] = sample.mass.Clone()
            
        h = samples[0].mass
        summc = ROOT.TH1F(cuts + dilepton + '_sumMC', '', h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        if sumMCfile is not None:
            summc.SetDirectory(sumMCfile)
        
        ## Sort by increasing integral.
        #samples.sort(key=lambda x: x.mass.Integral(x.mass.FindBin(150), x.mass.FindBin(1e9)))
        #print dilepton, [x.name for x in samples]

        hdata = data[dilepton]

        # Print a pretty table.
        for mass_range in [(60,120), (120,200), (120,), (200,)]:
            print 'cuts: %s  dilepton: %s  mass range: %s' % (cuts, dilepton, mass_range)
            for sample in samples:
                sample.integral = get_integral(sample.mass, *mass_range, integral_only=True, include_last_bin=False)
            hdata_integral = get_integral(hdata, *mass_range, integral_only=True, include_last_bin=False)
            print '%100s%20s%20s +/- %20s%20s' % ('sample', 'weight for %i/nb' % int(int_lumi*1000), 'integral', 'error', 'limit if int=0')
            print '%100s%20s%20.6f +/- %20.6f%20s' % ('data', '-', hdata_integral, hdata_integral**0.5, '-')
            sum_mc = 0.
            var_sum_mc = 0.
            sums = defaultdict(float)
            var_sums = defaultdict(float)
            for sample in sorted(samples, key=lambda x: x.integral, reverse=True):
                w = sample.partial_weight*int_lumi
                var = w * sample.integral # not w**2 * sample.integral because sample.integral is already I*w
                if 'zssm' not in sample.name:
                    sum_mc += sample.integral
                    var_sum_mc += var
                if joins.has_key(sample.name):
                    join_name = joins[sample.name]
                    sums[join_name] += sample.integral
                    var_sums[join_name] += var
                if sample.integral == 0:
                    limit = '%.6f' % (3*w)
                else:
                    limit = '-'
                print '%100s%20.6f%20.6f +/- %20.6f%20s' % (sample.nice_name, w, sample.integral, var**0.5, limit)
            print '%100s%20s%20.6f +/- %20.6f%20s' % ('sum MC (not including Z\')', '-', sum_mc, var_sum_mc**0.5, '-')
            for join_name in sorted(sums.keys()):
                print '%100s%20s%20.6f +/- %20.6f%20s' % (join_name, '-', sums[join_name], var_sums[join_name]**0.5, '-')
            print
        print
        
        # Stack 'em.
        s = ROOT.THStack('s', '')

        last_mc = None
        for sample in samples:
            h = sample.mass
            color = joins_colors[joins[sample.name]] if do_joins and joins.has_key(sample.name) else sample.color
            h.SetLineColor(color)
            h.SetMarkerStyle(0)
            if 'zssm' not in sample.name:
                h.SetFillColor(color)
                s.Add(h)
                summc.Add(h)
            last_mc = h

        if draw_zssm and dilepton == 'MuonsPlusMuonsMinus':
            l = ROOT.TLegend(0.63, 0.54, 0.87, 0.87)
        else:
            l = ROOT.TLegend(0.70, 0.59, 0.87, 0.87)
        l.SetFillColor(0)

        m = ROOT.TMarker()
        m.SetMarkerStyle(20)
        m.SetMarkerSize(0.5)
        m.SetMarkerColor(ROOT.kBlack)
        l.AddEntry(m, 'Data', 'LPE')

        legend_already = set()
        for sample in reversed(samples):
            nice_name = sample.nice_name
            if do_joins and joins.has_key(sample.name):
                join_nice_name = joins[sample.name]
                if join_nice_name in legend_already:
                    continue
                else:
                    legend_already.add(join_nice_name)
                    nice_name = join_nice_name
            if 'zssm' in sample.name and (not draw_zssm or dilepton != 'MuonsPlusMuonsMinus'):
                continue
            l.AddEntry(sample.mass, nice_name, 'F')

        s.Draw('hist')
        # must call Draw first or the THStack doesn't have a histogram/axis
        s.GetXaxis().SetRangeUser(*xax)

        hdata = data[dilepton]
        hdata.Rebin(rebin_factor)
        move_above_into_bin(hdata, xax[1])

        mymin = real_hist_min(s.GetStack().Last(), user_range=xax) * 0.7
        #mymin = real_hist_min(last_mc, user_range=xax) * 0.7
        mymax = real_hist_max(s.GetStack().Last(), user_range=xax, use_error_bars=False) * 1.05
        if hdata.GetEntries() > 0:
            rhm = real_hist_max(hdata, user_range=xax)
            mymax = max(mymax, rhm)

        if use_yaxis and yaxis.has_key(dilepton):
            yaxismin, yaxismax = yaxis[dilepton]
            if yaxismin is not None:
                mymin = yaxismin
            if yaxismax is not None:
                mymax = yaxismax

        #sys.stderr.write('%s %s (real s min %s) %s %s\n' % ( cuts, dilepton, s.GetMinimum(), mymin, mymax))
        
        s.SetMinimum(mymin)
        s.SetMaximum(mymax)

        hdata.SetTitle('')
        hdata.GetXaxis().SetRangeUser(*xax)
        hdata.GetXaxis().SetTitle(titleize[to_compare] % (subtitleize[dilepton], unitize[to_compare]))
        hdata.GetYaxis().SetTitle('Events/%i %s' % (rebin_factor, unitize[to_compare]))
        hdata.GetYaxis().SetTitleOffset(1.2)
        hdata.SetMinimum(mymin)
        hdata.SetMaximum(mymax)
        hdata.SetMarkerStyle(20)
        hdata.SetMarkerSize(0.5)
        hdata.SetStats(0)
        hdata.Draw('same e1')

        if draw_zssm and dilepton == 'MuonsPlusMuonsMinus':
            from samples import zssm750
            zp = zssm750.mass
            zp.SetTitle('')
            zp.SetLineWidth(2)
            zp.GetXaxis().SetRangeUser(*xax)
            zp.SetMinimum(mymin)
            zp.SetMaximum(mymax)
            zp.SetStats(0)
            zp.Draw('hist same')
        
        t1 = ROOT.TLatex(0.4, 0.93, '#sqrt{s} = 7 TeV,  #int L dt = %.1f pb^{-1}' % int_lumi)
        t2 = ROOT.TLatex(0.1, 0.93, 'CMS preliminary')

        for t in t1, t2:
            t.SetTextSize(0.0375)
            t.SetNDC()

        t1.Draw()
        t2.Draw()

        l.Draw('same')

        data_for_file = hdata.Clone(cuts + dilepton + '_data')
        if sumMCfile is not None:
            data_for_file.SetDirectory(sumMCfile)
        
        ps.save(dilepton)

        if False and dilepton == 'MuonsPlusMuonsMinus':
            ps.c.SetLogx(1)
            ps.save(dilepton + '_logx')
            ps.c.SetLogx(0)

        data_c = cumulative_histogram(hdata)
        summc_c = cumulative_histogram(summc)
        summc_c.SetLineColor(4)
        data_c.SetMarkerStyle(20)
        data_c.SetMarkerSize(0.5)

        for h in (data_c, summc_c):
            h.SetLineWidth(2)
            h.GetXaxis().SetRangeUser(*x_axis_limits2)
            h.SetStats(0)
            h.SetTitle('')
            #h.GetYaxis().SetRangeUser(0.1, 6.5e3)
            h.GetXaxis().SetTitle(titleize[to_compare] % (subtitleize[dilepton], unitize[to_compare]))
            h.GetYaxis().SetTitle('Events #geq X')
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetYaxis().SetLabelSize(0.028)

        if data_c.GetMaximum() > summc_c.GetMaximum():
            data_c.Draw('pe')
            summc_c.Draw('same')
        else:
            summc_c.Draw()
            data_c.Draw('pe same')

        t1.Draw()
        t2.Draw()

        l = ROOT.TLegend(0.60, 0.69, 0.86, 0.85)
        l.SetFillColor(0)
        l.AddEntry(data_c, 'Data', 'LPE')
        l.AddEntry(summc_c, 'Bkg. prediction', 'L')
        l.Draw('same')

        ps.save(dilepton + '_cumulative')

        if sumMCfile is not None:
            sumMCfile.Write()

    if False:
        for sample in samples:
            h = sample.histos[cuts + 'MuonsPlusMuonsMinus'].Clone()
            h.Scale(1/sample.partial_weight/int_lumi)
            h2 = sample.histos[cuts + 'MuonsElectronsOppSign'].Clone()
            h2.Scale(1/sample.partial_weight/int_lumi)
            h.Divide(h2)
            h.Draw('hist')
            ps.save('mumuemuratio_%s' % sample.name)

if hadd_tmp:
    os.system('rm %s' % data_fn)
