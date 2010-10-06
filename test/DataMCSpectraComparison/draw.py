#!/usr/bin/env python

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import real_hist_max, real_hist_min, set_zp2mu_style, plot_saver, ROOT
set_zp2mu_style()

# 166.9 nb in jul15.root, 2880.2 in prompt.root
int_lumi = (166.9 + 2880.2)/1000 # /pb
#int_lumi = 100

from samples import *

rebin_factor = 10
x_axis_limits = 110, 400
to_compare = 'DileptonMass'

for x in sys.argv:
    if 'xax' in x:
        xax0, xax1 = x.split('=')[1].split(',')
        x_axis_limits = int(xax0), int(xax1)
    if 'compare' in x:
        to_compare = x.split('=')[1]

histo_dir = 'ana_datamc'
def dir_name(c, d):
    return c + d + 'Histos'

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
    'DileptonMass': 'M_{%s} (%s)',
    'DileptonPt': '%s p_{T} (%s)'
    }
unitize = {
    'DileptonMass': 'GeV',
    'DileptonPt': 'GeV'
    }

#dileptons = ['MuonsPlusMuonsMinus', 'MuonsPlusMuonsPlus', 'MuonsMinusMuonsMinus', 'MuonsSameSign', 'ElectronsPlusElectronsMinus', 'ElectronsPlusElectronsPlus', 'ElectronsMinusElectronsMinus', 'ElectronsSameSign', 'MuonsPlusElectronsMinus', 'MuonsMinusElectronsPlus', 'MuonsPlusElectronsPlus', 'MuonsMinusElectronsMinus', 'MuonsElectronsOppSign', 'MuonsElectronsSameSign']
#cutss = ['Std','VBTF','Pt20']

dileptons = ['MuonsPlusMuonsMinus', 'MuonsSameSign', 'ElectronsPlusElectronsMinus', 'ElectronsSameSign', 'MuonsElectronsOppSign']
cutss = ['VBTF']

ROOT.TH1.AddDirectory(False)

ps = plot_saver('plots/datamc')

def integ(h,a,b=1e9):
    return h.Integral(h.FindBin(a), h.FindBin(b))

#samples = [s for s in samples if not s.name in ['ww', 'zz', 'wz', 'qcd500']]

for cuts in cutss:
    plot_dir = 'plots/datamc/%s/%i_%i/%s' % (to_compare, x_axis_limits[0], x_axis_limits[1], cuts)
    ps.set_plot_dir(plot_dir)

    fdata = ROOT.TFile(os.path.join(histo_dir, 'ana_datamc_data.root'))
    data = dict((d, getattr(fdata, dir_name(cuts, d)).Get(to_compare).Clone()) for d in dileptons)

    for dilepton in dileptons:
        for sample in samples:
            # It would be more efficient to have the sample loop be
            # the outer one thanks to the file opening/closing, but
            # the code is cleaner this way.
            f = ROOT.TFile(os.path.join(histo_dir, 'ana_datamc_%s.root' % sample.name))
            sample.mass = getattr(f, dir_name(cuts, dilepton)).Get(to_compare).Clone()
            sample.mass.Rebin(rebin_factor)
            sample.mass.Scale(sample.partial_weight * int_lumi)

        ## Sort by increasing integral.
        #samples.sort(key=lambda x: x.mass.Integral(x.mass.FindBin(150), x.mass.FindBin(1e9)))
        #print dilepton, [x.name for x in samples]

        hdata = data[dilepton]

        for mass_range in [(60,120), (120,)]:
            print 'cuts: %s  dilepton: %s  mass range: %s' % (cuts, dilepton, mass_range)
            for sample in samples:
                sample.integral = integ(sample.mass, *mass_range)
            hdata_integral = integ(hdata, *mass_range)
            print '%100s%20s%20s +/- %20s%20s' % ('sample', 'weight for %i/nb' % int(int_lumi*1000), 'integral', 'error', 'limit if int=0')
            print '%100s%20s%20.6f +/- %20.6f%20s' % ('data', '-', hdata_integral, hdata_integral**0.5, '-')
            sum_mc = 0.
            var_sum_mc = 0.
            for sample in sorted(samples, key=lambda x: x.integral, reverse=True):
                w = sample.partial_weight*int_lumi
                sum_mc += sample.integral
                var_sum_mc += w**2 * sample.integral
                if sample.integral == 0:
                    limit = '%.6f' % (3*w)
                else:
                    limit = '-'
                print '%100s%20.6f%20.6f +/- %20.6f%20s' % (sample.nice_name, w, sample.integral, w*sample.integral**0.5, limit)
            print '%100s%20s%20.6f +/- %20.6f%20s' % ('sum MC', '-', sum_mc, var_sum_mc**0.5, '-')
            print
        print
        
        # Stack 'em.
        s = ROOT.THStack('s', '')

        last_mc = None
        for sample in samples:
            if 'qcd' in sample.name:
                continue
            h = sample.mass
            h.SetFillColor(sample.color)
            h.SetLineColor(sample.color)
            h.SetMarkerStyle(0)
            s.Add(h)
            last_mc = h

        l = ROOT.TLegend(0.69, 0.56, 0.87, 0.87)
        l.SetFillColor(0)

        m = ROOT.TMarker()
        m.SetMarkerStyle(20)
        m.SetMarkerSize(0.5)
        m.SetMarkerColor(ROOT.kBlack)
        l.AddEntry(m, 'Data')

        for sample in reversed(samples):
            if 'qcd' in sample.name:
                continue
            l.AddEntry(sample.mass, sample.nice_name, 'F')

        s.Draw('hist')
        # must call Draw first or the THStack doesn't have a histogram/axis
        s.GetXaxis().SetRangeUser(*x_axis_limits)

        hdata = data[dilepton]
        hdata.Rebin(rebin_factor)

        mymin = real_hist_min(s.GetStack().Last(), user_range=x_axis_limits) * 0.7
        #mymin = real_hist_min(last_mc, user_range=x_axis_limits) * 0.7
        mymax = real_hist_max(s.GetStack().Last(), user_range=x_axis_limits, use_error_bars=False) * 1.05
        if hdata.GetEntries() > 0:
            rhm = real_hist_max(hdata, user_range=x_axis_limits)
            mymax = max(mymax, rhm)

        #sys.stderr.write('%s %s (real s min %s) %s %s\n' % ( cuts, dilepton, s.GetMinimum(), mymin, mymax))
        
        s.SetMinimum(mymin)
        s.SetMaximum(mymax)
        
        hdata.SetTitle('')
        hdata.GetXaxis().SetRangeUser(*x_axis_limits)
        hdata.GetXaxis().SetTitle(titleize[to_compare] % (subtitleize[dilepton], unitize[to_compare]))
        hdata.GetYaxis().SetTitle('Events/10 %s' % unitize[to_compare])
        hdata.SetMinimum(mymin)
        hdata.SetMaximum(mymax)
        hdata.SetMarkerStyle(20)
        hdata.SetMarkerSize(0.5)
        hdata.SetStats(0)
        hdata.Draw('same e1')

        t1 = ROOT.TLatex(0.4, 0.93, '#sqrt{s} = 7 TeV,  #int L dt = %.1f pb^{-1}' % int_lumi)
        t2 = ROOT.TLatex(0.1, 0.93, 'CMS preliminary')

        for t in t1, t2:
            t.SetTextSize(0.0375)
            t.SetNDC()

        t1.Draw()
        t2.Draw()

        l.Draw('same')

        ps.save(dilepton)

#        raise 1
