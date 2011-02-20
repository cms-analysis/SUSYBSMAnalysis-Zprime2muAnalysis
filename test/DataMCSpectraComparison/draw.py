#!/usr/bin/env python

# (py draw.py ana_datamc_current/muonsonly >! out.draw.muonsonly) && (py draw.py ana_datamc_current/allgood >! out.draw.allgood) && mv out.draw.* plots/datamc_current/ && tlock ~/asdf/plots.tgz plots/datamc_current

import sys, os, glob
from collections import defaultdict
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import cumulative_histogram, get_integral, move_above_into_bin, plot_saver, poisson_intervalize, real_hist_max, real_hist_min, set_zp2mu_style, ROOT

set_zp2mu_style()
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadRightMargin(0.07)

from samples import *

rebin_factor = 5
x_axis_limits = 50, 1050
x_axis_limits2 = 50, 500

to_compare = 'DileptonMass'
global_rescale = 3273/3404.6 if False else None
draw_zssm = True
use_poisson_intervals = True

do_joins = True
joins = [(s.name, 'jets') for s in samples if 'qcd' in s.name]
joins += [(x, 'jets') for x in ['inclmu15', 'wmunu', 'wjets']]
#joins += [(x, 'other prompt leptons') for x in ['singletop_tW', 'ztautau', 'ww', 'wz', 'zz']]
joins += [(x, 't#bar{t} + other prompt leptons') for x in ['ttbar', 'singletop_tW', 'ztautau', 'ww', 'wz', 'zz']]
joins += [(s.name, '#gamma/Z #rightarrow #mu^{+}#mu^{-}') for s in samples if 'dy' in s.name]
joins += [('zmumu', '#gamma/Z #rightarrow #mu^{+}#mu^{-}')]
joins = dict(joins)
joins_colors = {'jets': 4, 't#bar{t} + other prompt leptons': 2, '#gamma/Z #rightarrow #mu^{+}#mu^{-}': 7}

histo_dir = [x for x in sys.argv if os.path.isdir(x)][0]

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
lumi_syst_frac = 0.11

subtitleize = {
    'MuonsPlusMuonsMinus': '#mu^{+}#mu^{-}',
    'MuonsPlusMuonsPlus':  '#mu^{+}#mu^{+}',
    'MuonsMinusMuonsMinus': '#mu^{-}#mu^{-}',
    'MuonsSameSign': '#mu^{#pm}#mu^{#pm}',
    'MuonsAllSigns': '#mu#mu',
    'ElectronsPlusElectronsMinus': 'e^{+}e^{-}',
    'ElectronsPlusElectronsPlus': 'e^{+}e^{+}',
    'ElectronsMinusElectronsMinus': 'e^{-}e^{-}',
    'ElectronsSameSign': 'e^{#pm}e^{#pm}',
    'ElectronsAllSigns': 'ee',
    'MuonsPlusElectronsMinus': '#mu^{+}e^{-}',
    'MuonsMinusElectronsPlus': '#mu^{-}e^{+}',
    'MuonsPlusElectronsPlus': '#mu^{+}e^{+}',
    'MuonsMinusElectronsMinus': '#mu^{-}e^{-}',
    'MuonsElectronsOppSign': '#mu^{+}e^{-}/#mu^{-}e^{+}',
    'MuonsElectronsSameSign': '#mu^{#pm}e^{#pm}',
    'MuonsElectronsAllSigns': 'e#mu',
    }
titleize = {
    'DileptonMass': 'm(%s)%s',
    'DileptonPt': '%s p_{T}%s',
    }
unitize = {
    'DileptonMass': ' (GeV)',
    'DileptonPt': ' (GeV)',
    }
yaxis = {
    ('MuonsPlusMuonsMinus', False): (4e-3, None),
    ('MuonsPlusMuonsMinus', True): (1.6, None),
#   ('MuonsSameSign', False): (5e-5, 2.5),
   ('MuonsElectronsOppSign', False): (8e-3, 30),
    }
use_yaxis = True

dileptons = ['MuonsPlusMuonsMinus', 'MuonsSameSign', 'MuonsAllSigns', 'MuonsElectronsOppSign', 'MuonsElectronsSameSign', 'MuonsElectronsAllSigns']
cutss = ['VBTF', 'Our', 'OurNoIso']

ROOT.TH1.AddDirectory(False)

def dir_name(c, d):
    return c + d + 'Histos'

pdir = 'plots/datamc'
if histo_dir != 'ana_datamc':
    pdir += '_' + histo_dir.replace('ana_datamc_', '')
ps = plot_saver(pdir, size=(600,600), pdf_log=True)
save_plots = 'no_plots' not in sys.argv

if global_rescale is not None:
    for s in samples:
        s.partial_weight *= global_rescale

#samples = [s for s in samples if not s.name in ['ww', 'zz', 'wz', 'qcd500']]

for cuts in cutss:
    plot_dir = pdir + '/%s/%s' % (to_compare, cuts)
    ps.set_plot_dir(plot_dir)
    for cumulative in (False, True):
        data = dict((d, getattr(fdata, dir_name(cuts, d)).Get(to_compare).Clone()) for d in dileptons)

        for dilepton in dileptons:
            if int_lumi > 39 and 'Electron' in dilepton:
                continue

            xax = x_axis_limits if (dilepton == 'MuonsPlusMuonsMinus' and not cumulative) else x_axis_limits2

            for sample in samples:
                # It would be more efficient to have the sample loop be
                # the outer one thanks to the file opening/closing, but
                # the code is cleaner this way.
                f = ROOT.TFile(os.path.join(histo_dir, 'ana_datamc_%s.root' % sample.name))
                sample.mass = getattr(f, dir_name(cuts, dilepton)).Get(to_compare).Clone()
                sample.mass.Rebin(rebin_factor)
                sample.mass.Scale(sample.partial_weight * int_lumi)
                if cumulative:
                    sample.mass_not_cumulative = sample.mass
                    sample.mass = cumulative_histogram(sample.mass)
                else:
                    move_above_into_bin(sample.mass, xax[1])

            ## Sort by increasing integral.
            #if not cumulative:
            #    samples.sort(key=lambda x: x.mass.Integral(x.mass.FindBin(150), x.mass.FindBin(1e9)))
            #    print dilepton, [x.name for x in samples]

            hdata = data[dilepton]

            # Print a pretty table.
            if not cumulative:
                for mass_range in [(60,120), (120,200), (120,), (200,), (586, 914), (150,200)]:
                    if mass_range == (586,914) and (cuts != 'Our' or dilepton != 'MuonsPlusMuonsMinus'):
                        continue
                    print 'cuts: %s  dilepton: %s  mass range: %s' % (cuts, dilepton, mass_range)
                    for sample in samples:
                        sample.integral = get_integral(sample.mass, *mass_range, integral_only=True, include_last_bin=False)
                    hdata_integral = get_integral(hdata, *mass_range, integral_only=True, include_last_bin=False)
                    print '%50s%20s%20s%20s%20s%20s%20s%20s' % ('sample', 'weight for %i/nb' % int(int_lumi*1000), 'integral', 'stat error', 'limit if int=0', 'syst error', 'lumi error', 'total error')
                    print '%50s%20s%20.6f%20.6f' % ('data', '-', hdata_integral, hdata_integral**0.5)
                    sum_mc = 0.
                    var_sum_mc = 0.
                    syst_var_sum_mc = 0.
                    sums = defaultdict(float)
                    var_sums = defaultdict(float)
                    syst_var_sums = defaultdict(float)
                    for sample in sorted(samples, key=lambda x: x.integral, reverse=True):
                        w = sample.partial_weight*int_lumi
                        var = w * sample.integral # not w**2 * sample.integral because sample.integral is already I*w
                        syst_var = (sample.syst_frac * sample.integral)**2

                        if 'zssm' not in sample.name:
                            sum_mc += sample.integral
                            var_sum_mc += var
                            syst_var_sum_mc += syst_var

                        if joins.has_key(sample.name):
                            join_name = joins[sample.name]
                            sums[join_name] += sample.integral
                            var_sums[join_name] += var
                            syst_var_sums[join_name] += syst_var

                        limit = '%.6f' % (3*w) if sample.integral == 0 else '-'

                        lumi_err = lumi_syst_frac * sample.integral
                        tot_err = (var + syst_var + lumi_err**2)**0.5

                        print '%50s%20.6f%20.6f%20.6f%20s%20.6f%20.6f%20.6f' % (sample.nice_name, w, sample.integral, var**0.5, limit, syst_var**0.5, lumi_err, tot_err)

                    lumi_err = lumi_syst_frac * sum_mc
                    tot_err = (var_sum_mc + syst_var_sum_mc + lumi_err**2)**0.5
                    print '%50s%20s%20.6f%20.6f%20s%20.6f%20.6f%20.6f' % ('sum MC (not including Z\')', '-', sum_mc, var_sum_mc**0.5, '-', syst_var_sum_mc**0.5, lumi_err, tot_err)
                    for join_name in sorted(sums.keys()):
                        lumi_err = lumi_syst_frac * sums[join_name]
                        tot_err = (var_sums[join_name] + syst_var_sums[join_name] + lumi_err**2)**0.5
                        print '%50s%20s%20.6f%20.6f%20s%20.6f%20.6f%20.6f' % (join_name, '-', sums[join_name], var_sums[join_name]**0.5, '-', syst_var_sums[join_name]**0.5, lumi_err, tot_err)
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
                if 'zssm' not in sample.name and ('MuonsElectrons' not in dilepton or ('zmumu' not in sample.name and 'dy' not in sample.name)):
                    h.SetFillColor(color)
                    s.Add(h)
                last_mc = h

            if draw_zssm and (dilepton == 'MuonsPlusMuonsMinus' and not cumulative):
                l = ROOT.TLegend(0.47, 0.55, 0.88, 0.88)
            elif dilepton == 'MuonsPlusMuonsMinus' and cumulative:
                l = ROOT.TLegend(0.47, 0.55, 0.88, 0.88)
            else:
                l = ROOT.TLegend(0.53, 0.69, 0.76, 0.84)
            l.SetFillColor(0)
            l.SetBorderSize(0)

            m = ROOT.TMarker()
            m.SetMarkerStyle(20)
            m.SetMarkerSize(1.0)
            m.SetMarkerColor(ROOT.kBlack)
            l.AddEntry(m, 'DATA', 'LP')

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
                if 'zssm' in sample.name and (not draw_zssm or dilepton != 'MuonsPlusMuonsMinus' or cumulative) or ('MuonsElectrons' in dilepton and ('zmumu' in sample.name or 'dy' in sample.name)):
                    continue
                l.AddEntry(sample.mass, nice_name, 'F')

            s.Draw('hist')
            # must call Draw first or the THStack doesn't have a histogram/axis
            s.GetXaxis().SetRangeUser(*xax)
            
            hdata.Rebin(rebin_factor)
            if cumulative:
                hdata = cumulative_histogram(hdata)
            else:
                move_above_into_bin(hdata, xax[1])

            mymin = real_hist_min(s.GetStack().Last(), user_range=xax) * 0.7
            #mymin = real_hist_min(last_mc, user_range=xax) * 0.7
            mymax = real_hist_max(s.GetStack().Last(), user_range=xax, use_error_bars=False) * 1.05
            if hdata.GetEntries() > 0:
                rhm = real_hist_max(hdata, user_range=xax)
                mymax = max(mymax, rhm)

            if use_yaxis and yaxis.has_key((dilepton,cumulative)):
                yaxismin, yaxismax = yaxis[(dilepton,cumulative)]
                if yaxismin is not None:
                    mymin = yaxismin
                if yaxismax is not None:
                    mymax = yaxismax

            s.SetMinimum(mymin)
            s.SetMaximum(mymax)

            hdata.SetStats(0)
            data_draw_cmd = 'same p e'
            if use_poisson_intervals:
                hdata = poisson_intervalize(hdata, True)
                data_draw_cmd += ' z'
            hdata.SetTitle('')
            hdata.GetXaxis().SetRangeUser(*xax)
            hdata.GetXaxis().SetTitle(titleize[to_compare] % (subtitleize[dilepton], unitize[to_compare]))
            if cumulative:
                title = 'Events #geq %s' % (titleize[to_compare] % (subtitleize[dilepton], ''))
            else:
                title = 'Events / %i%s' % (rebin_factor, unitize[to_compare].replace('(','').replace(')',''))
            hdata.GetYaxis().SetTitle(title)
            hdata.GetXaxis().SetTitleOffset(0.9)
            hdata.GetXaxis().SetTitleSize(0.047)
            hdata.GetYaxis().SetTitleOffset(1.2)
            hdata.GetYaxis().SetTitleSize(0.047)
            hdata.SetMinimum(mymin)
            hdata.SetMaximum(mymax)
            hdata.SetMarkerStyle(20)
            hdata.SetMarkerSize(1.0)
            hdata.Draw(data_draw_cmd)

            if draw_zssm and not cumulative and dilepton == 'MuonsPlusMuonsMinus':
                from samples import zssm750
                zp = zssm750.mass
                zp.SetTitle('')
                zp.SetLineWidth(2)
                zp.GetXaxis().SetRangeUser(*xax)
                zp.SetMinimum(mymin)
                zp.SetMaximum(mymax)
                zp.SetStats(0)
                zp.Draw('hist same')

            t1 = ROOT.TPaveLabel(0.458, 0.893, 0.908, 0.993, '#sqrt{s} = 7 TeV, #int L dt = %.f pb^{-1}' % round(int_lumi), 'brNDC')
            t2 = ROOT.TPaveLabel(0.225, 0.907, 0.332, 0.997, 'CMS preliminary', 'brNDC')
            t1.SetTextSize(0.35)
            t2.SetTextSize(0.45)
            for t in t1, t2:
                t.SetBorderSize(0)
                t.SetFillColor(0)
                t.SetFillStyle(0)
                t.Draw()

            l.SetTextSize(0.03)
            l.Draw('same')

            n = dilepton
            if cumulative:
                n += '_cumulative'
            if save_plots:
                ps.save(n)

if hadd_tmp:
    os.system('rm %s' % data_fn)
