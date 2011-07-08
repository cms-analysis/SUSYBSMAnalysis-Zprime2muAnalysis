#!/usr/bin/env python

# (py draw.py ana_datamc_Run2011AMuonsOnly >! plots/out.draw.ana_datamc_Run2011AMuonsOnly) && (py draw.py ana_datamc_Run2011AMuonsOnly compare=DimuonMassVertexConstrained >! plots/out.draw.ana_datamc_Run2011AMuonsOnly_DimuonMassVertexConstrained) && (py draw.py ana_datamc_Run2011A >! plots/out.draw.ana_datamc_Run2011A) && tlp plots/datamc* plots/out.draw.*

import sys, os, glob
from pprint import pprint
from collections import defaultdict
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import cumulative_histogram, get_integral, move_above_into_bin, plot_saver, poisson_intervalize, real_hist_max, real_hist_min, set_zp2mu_style, ROOT

set_zp2mu_style()
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadRightMargin(0.07)

from samples import *

rebin_factor = 5
x_axis_limits = 70, 1300
x_axis_limits2 = 70, 1300

to_compare = 'DileptonMass'
for x in sys.argv:
    if x.startswith('compare='):
        to_compare = x.replace('compare=', '')
        break

global_rescale = {
    'VBTF': 236163/202178.739663,
    'OurNew': 275060/237045.128578,
    'OurOld': 284999/245756.305310,
    'OurNoIso': 279681/240450.139476,
    }
if 'norescale' in sys.argv:
    global_rescale = {}
def get_rescale_factor(cuts, dilepton):
    # Don't rescale e-mu plots.
    return 1. if 'electron' in dilepton.lower() else global_rescale.get(cuts, 1.)
    
draw_zssm = True 
use_poisson_intervals = True
overflow_bin = True

do_joins = True
joins = [(s.name, 'jets') for s in samples if 'qcd' in s.name]
joins += [(x, 'jets') for x in ['inclmu15', 'wmunu', 'wjets']]
#joins += [(x, 'other prompt leptons') for x in ['singletop_tW', 'ztautau', 'ww', 'wz', 'zz']]
joins += [(x, 't#bar{t} + other prompt leptons') for x in ['ttbar', 'singletop_tW', 'ztautau', 'ww', 'wz', 'zz']]
joins += [(s.name, '#gamma/Z #rightarrow #mu^{+}#mu^{-}') for s in samples if 'dy' in s.name]
joins += [('zmumu', '#gamma/Z #rightarrow #mu^{+}#mu^{-}')]
joins = dict(joins)
joins_colors = {'jets': 4, 't#bar{t} + other prompt leptons': 2, 'other prompt leptons': 2, '#gamma/Z #rightarrow #mu^{+}#mu^{-}': 7}

histo_dir = [x for x in sys.argv if os.path.isdir(x)][0]

data_fn = os.path.join(histo_dir, 'ana_datamc_data.root')
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

int_lumi = parse_lumi_from_log(data_fn.replace('.root', '.lumi'))
lumi_syst_frac = 0.06

print '"joins" are:'
pprint(joins)
print 'total lumi from data: %.1f/pb' % int_lumi
print 'rescaling mumu MC histograms by these cut-dependent factors:'
pprint(global_rescale)
print 'comparing', to_compare
print 'using poisson error bars on plots:', use_poisson_intervals
print 'last bin contains overflow:', overflow_bin

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
    'DimuonMassVertexConstrained': 'm(%s)%s',
    'DileptonPt': '%s p_{T}%s',
    }
unitize = {
    'DileptonMass': ' [GeV]',
    'DimuonMassVertexConstrained': ' [GeV]',
    'DileptonPt': ' [GeV]',
    }
yaxis = {
    #('MuonsSameSign', False): (1e-4, 5),
    }
use_yaxis = True

dileptons = ['MuonsPlusMuonsMinus', 'MuonsSameSign', 'MuonsAllSigns', 'MuonsElectronsOppSign', 'MuonsElectronsSameSign', 'MuonsElectronsAllSigns']
cutss = ['VBTF', 'OurNew', 'OurOld', 'Simple', 'EmuVeto', 'OurNoIso']
#cutss = ['OurNew']
mass_ranges_for_table = [(60,120), (120,200), (200,), (200,400), (400,), (600,)]

if 'forscale' in sys.argv:
    global_rescale = {}
    mass_ranges_for_table = [(60,120)]
    dileptons = ['MuonsPlusMuonsMinus']
    cutss.remove('Simple')
    cutss.remove('EmuVeto')

ROOT.TH1.AddDirectory(False)

def dir_name(c, d):
    return c + d + 'Histos'

pdir = 'plots/datamc'
if histo_dir != 'ana_datamc':
    pdir += '_' + histo_dir.replace('ana_datamc_', '')
ps = plot_saver(pdir, size=(900,600), pdf_log=True, pdf=True)
save_plots = 'no_plots' not in sys.argv

#samples = [s for s in samples if not s.name in ['wjets']]

for cuts in cutss:
    if not hasattr(fdata, dir_name(cuts, 'MuonsPlusMuonsMinus')):
        continue
    plot_dir = pdir + '/%s/%s' % (to_compare, cuts)
    ps.set_plot_dir(plot_dir)

    if cuts == 'EmuVeto':
        dils = [x for x in dileptons if 'Electron' in x]
    else:
        dils = dileptons

    if 'Dimuon' in to_compare:
        dils = [x for x in dils if 'Electron' not in x]

    for cumulative in (False, True):
        data = dict((d, getattr(fdata, dir_name(cuts, d)).Get(to_compare).Clone()) for d in dils)

        for dilepton in dils:
            xax = x_axis_limits if (dilepton == 'MuonsPlusMuonsMinus' and not cumulative) else x_axis_limits2

            for sample in samples:
                # It would be more efficient to have the sample loop be
                # the outer one thanks to the file opening/closing, but
                # the code is cleaner this way.
                f = ROOT.TFile(os.path.join(histo_dir, 'mc', 'ana_datamc_%s.root' % sample.name))
                sample.mass = getattr(f, dir_name(cuts, dilepton)).Get(to_compare).Clone()
                sample.mass.Rebin(rebin_factor)
                sample.scaled_by = get_rescale_factor(cuts, dilepton) * sample.partial_weight * int_lumi
                sample.mass.Scale(sample.scaled_by)
                if cumulative:
                    sample.mass_not_cumulative = sample.mass
                    sample.mass = cumulative_histogram(sample.mass)
                elif overflow_bin:
                    move_above_into_bin(sample.mass, xax[1])

            ## Sort by increasing integral.
            #if not cumulative:
            #    samples.sort(key=lambda x: x.mass.Integral(x.mass.FindBin(150), x.mass.FindBin(1e9)))
            #    print dilepton, [x.name for x in samples]

            hdata = data[dilepton]

            # Print a pretty table.
            if not cumulative:
                for mass_range in mass_ranges_for_table:
                    if mass_range == (586,914) and (cuts != 'Our' or dilepton != 'MuonsPlusMuonsMinus'):
                        continue
                    print 'cuts: %s  dilepton: %s  mass range: %s' % (cuts, dilepton, mass_range)
                    for sample in samples:
                        sample.integral = get_integral(sample.mass, *mass_range, integral_only=True, include_last_bin=False)
                        sample.raw_integral = sample.integral / sample.scaled_by
                    hdata_integral = get_integral(hdata, *mass_range, integral_only=True, include_last_bin=False)
                    add_lumi_error_to_total = 'Electron' in dilepton or not global_rescale.has_key(cuts)
                    print '%50s%20s%20s%20s%20s%20s%20s%20s%20s' % ('sample', 'weight for %.1f/pb' % int_lumi, 'raw integral', 'integral', 'stat error', 'limit if int=0', 'syst error', 'lumi error', 'tot error' + (' (nolumi)' if not add_lumi_error_to_total else ''))
                    print '%50s%20s%20i%20.6f%20.6f' % ('data', '-', int(hdata_integral), hdata_integral, hdata_integral**0.5)
                    sum_mc = 0.
                    var_sum_mc = 0.
                    syst_var_sum_mc = 0.
                    sums = defaultdict(float)
                    var_sums = defaultdict(float)
                    syst_var_sums = defaultdict(float)
                    for sample in sorted(samples, key=lambda x: x.integral, reverse=True):
                        w = sample.scaled_by
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
                        if add_lumi_error_to_total:
                            tot_err = (var + syst_var + lumi_err**2)**0.5
                        else:
                            tot_err = (var + syst_var)**0.5

                        print '%50s%20.6f%20f%20.6f%20.6f%20s%20.6f%20.6f%20.6f' % (sample.nice_name, w, sample.raw_integral, sample.integral, var**0.5, limit, syst_var**0.5, lumi_err, tot_err)

                    lumi_err = lumi_syst_frac * sum_mc
                    if add_lumi_error_to_total:
                        tot_err = (var_sum_mc + syst_var_sum_mc + lumi_err**2)**0.5
                    else:
                        tot_err = (var_sum_mc + syst_var_sum_mc)**0.5
                        
                    print '%50s%20s%20s%20.6f%20.6f%20s%20.6f%20.6f%20.6f' % ('sum MC (not including Z\')', '-', '-', sum_mc, var_sum_mc**0.5, '-', syst_var_sum_mc**0.5, lumi_err, tot_err)
                    for join_name in sorted(sums.keys()):
                        lumi_err = lumi_syst_frac * sums[join_name]
                        if add_lumi_error_to_total:
                            tot_err = (var_sums[join_name] + syst_var_sums[join_name] + lumi_err**2)**0.5
                        else:
                            tot_err = (var_sums[join_name] + syst_var_sums[join_name])**0.5
                        print '%50s%20s%20s%20.6f%20.6f%20s%20.6f%20.6f%20.6f' % (join_name, '-', '-', sums[join_name], var_sums[join_name]**0.5, '-', syst_var_sums[join_name]**0.5, lumi_err, tot_err)
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
                last_mc = h

            if draw_zssm and (dilepton == 'MuonsPlusMuonsMinus' and not cumulative):
                l = ROOT.TLegend(0.47, 0.55, 0.88, 0.88)
            elif dilepton == 'MuonsPlusMuonsMinus' and cumulative:
                l = ROOT.TLegend(0.47, 0.55, 0.88, 0.88)
            else:
                l = ROOT.TLegend(0.50, 0.69, 0.76, 0.88)
            l.SetFillColor(0)
            l.SetBorderSize(0)

            m = ROOT.TMarker()
            m.SetMarkerStyle(20)
            m.SetMarkerSize(0.8)
            m.SetMarkerColor(ROOT.kBlack)
            l.AddEntry(m, 'DATA', 'EP')

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

            xtitle = titleize[to_compare] % (subtitleize[dilepton], unitize[to_compare])
            if cumulative:
                ytitle = 'Events #geq %s' % (titleize[to_compare] % (subtitleize[dilepton], ''))
            else:
                ytitle = 'Events / %i%s' % (rebin_factor, unitize[to_compare].translate(None, '()[]'))

            s.SetTitle(';%s;%s' % (xtitle, ytitle))
            s.Draw('hist')
            # must call Draw first or the THStack doesn't have a histogram/axis
            s.GetXaxis().SetRangeUser(*xax)
            s.GetXaxis().SetTitleOffset(0.9)
            s.GetXaxis().SetTitleSize(0.047)
            s.GetYaxis().SetTitleOffset(1.2)
            s.GetYaxis().SetTitleSize(0.047)

            hdata.Rebin(rebin_factor)
            if cumulative:
                hdata = cumulative_histogram(hdata)
            elif overflow_bin:
                move_above_into_bin(hdata, xax[1])
            else:
                overflow_integral = get_integral(hdata, xax[1], integral_only=True)
                if overflow_integral > 0:
                    print 'WARNING: in %s, data histogram has points in overflow (mass bins above %.f GeV)! integral = %f' % (cuts + ' ' + dilepton,  xax[1], overflow_integral)

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
            hdata.GetXaxis().SetRangeUser(*xax)
            hdata.SetMinimum(mymin)
            hdata.SetMaximum(mymax)
            hdata.SetMarkerStyle(20)
            hdata.SetMarkerSize(0.8)
            hdata.Draw(data_draw_cmd)

            if draw_zssm and not cumulative and dilepton == 'MuonsPlusMuonsMinus':
                from samples import zssm1000
                zp = zssm1000.mass
                zp.SetTitle('')
                zp.SetLineWidth(2)
                zp.GetXaxis().SetRangeUser(*xax)
                zp.SetMinimum(mymin)
                zp.SetMaximum(mymax)
                zp.SetStats(0)
                zp.Draw('hist same')

            t = ROOT.TPaveLabel(0.20, 0.89, 0.86, 0.99, 'CMS 2011 preliminary   #sqrt{s} = 7 TeV    #int L dt = %.f pb^{-1}' % round(int_lumi), 'brNDC')
            t.SetTextSize(0.35)
            t.SetBorderSize(0)
            t.SetFillColor(0)
            t.SetFillStyle(0)
            t.Draw()

            l.SetTextSize(0.03)
            l.Draw('same')
            ## huge crappy hack for "EP" in TLegend::AddEntry not working
            #ll = ROOT.TLine()
            #ll.DrawLineNDC(0.532, 0.82, 0.532, 0.875)

            n = dilepton
            if cumulative:
                n += '_cumulative'
            if save_plots:
                ps.save(n)
