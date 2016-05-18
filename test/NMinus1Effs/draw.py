#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

# Import nm1entry.py
# - from MCSamples and roottools import *
# - nm1entry class
# - eff_wald()
# - nminus1s list
# - pretty dictionary
# - styles dictionary
from nm1entry import *

def draw_test(tag, printStats, lumi, do_tight=False):
    # Specific stylings for this script on top of zp2mu_style
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ROOT.gStyle.SetTitleX(0.25)
    ROOT.gStyle.SetTitleY(0.50)
    #ROOT.gStyle.SetTitleH(0.07)
    ROOT.TH1.AddDirectory(0)
    outfile = ROOT.TFile("whargl.root","recreate")
    iarp=0

    # 'plots' = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/NMinus1Effs/plots/TAG/tag'
    # TAG = MC/Data production TAG
    # tag = sub tag for each version of plots made with a single production TAG
    psn = 'plots'
    if tag:
        psn = 'plots/%s'%tag

    # cjsbad - fix the referencing to tightnm1...
    #if do_tight:
    #    psn += '_tight'
    #    nminus1s = tightnm1[:]

    ps = plot_saver(psn, size=(600,600), log=False, pdf=True)
    ps.c.SetBottomMargin(0.2)

    data = nm1entry('data', True, lumi)#lumiCD )

    # 50 < m < 120 GeV
    samples50m120 = [dy50to120,dy120to200,dy200to400,ttbar_pow,ww_incl,wz,zz_incl,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]
    mc_samples_50m120 = [nm1entry(sample,False,lumi) for sample in samples50m120]

    # m > 120 GeV
    samples120m = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,wz,zz_incl,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]
    mc_samples_120m = [nm1entry(sample,False,lumi) for sample in samples120m]

    mass_ranges = [
         ('60m120',(60,120)),
         ('120m',(120,2500)),
        ]

    to_use = {
        '60m120':   [mc_samples_50m120,data],
        '120m':   [mc_samples_120m,data],
        }

    ymin = {
        '60m120':  0.7,
        '120m':    0.5,
        }

    #global_ymin = 0.
    global_ymin = None


    for name, mass_range in mass_ranges:
        pretty_name = pretty[name]

        lg = ROOT.TLegend(0.25, 0.21, 0.91, 0.44)
        lg.SetTextSize(0.03)
        lg.SetFillColor(0)
        lg.SetBorderSize(1)
        
        same = 'A'
        effs = []
        
        for entry in to_use[name]:

            l = len(nminus1s)
            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                color, fill = styles[entry.name]
                nminus1_num = ROOT.TH1F('num', '', l, 0, l)
                nminus1_den = ROOT.TH1F('den', '', l, 0, l)
                hnum = entry.histos['NoNo']
                num = get_integral(hnum, *mass_range, integral_only=True, include_last_bin=False, nm1=True)
                for i,nminus1 in enumerate(nminus1s):
                    if printStats: table(entry,nminus1, mass_range)
                    hden = entry.histos[nminus1]
                    den = get_integral(hden, *mass_range, integral_only=True, include_last_bin=False, nm1=True)
                    nminus1_num.SetBinContent(i+1, num)
                    nminus1_den.SetBinContent(i+1, den)
                eff,ycp,eylcp,eyhcp = binomial_divide(nminus1_num, nminus1_den)
            else:
                color = ROOT.kGreen+2
                fill = 1001
                nMC = len(entry)
                num = []
                pw = []
                den = []
                for i, mc in enumerate(entry):
                    den.append([])
                    hnum = mc.histos['NoNo']
                    num.append(get_integral(hnum, *mass_range, integral_only=True, include_last_bin=False, nm1=True))
                    pw.append(mc.partial_weight * lumi)
                    for j, nminus1 in enumerate(nminus1s):
                        if printStats: table_wald(mc,nminus1, mass_range)
                        hden = mc.histos[nminus1]
                        den[i].append(get_integral(hden, *mass_range, integral_only=True, include_last_bin=False, nm1=True))
                nm1Temp = ROOT.TH1F('temp','',l,0,l)
                eff,y,eyl,eyh = eff_wald(num, den, l, nMC, pw, nm1Temp)
            eff.SetTitle(pretty_name)
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin[name], 1.01)
            eff.GetXaxis().SetTitle('cut')
            eff.GetYaxis().SetLabelSize(0.027)
            eff.GetYaxis().SetTitle('n-1 efficiency')
            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                draw = 'P'
                eff.SetLineColor(color)
                eff.SetMarkerStyle(20)
                eff.SetMarkerSize(1.05)
                eff.SetMarkerColor(color)
                #lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
                lg.AddEntry(eff, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')
            else:
                draw = '2'
                eff.SetLineColor(color)
                eff.SetFillColor(color)
                eff.SetFillStyle(fill)
                #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
                lg.AddEntry(eff, 'Simulation', 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            bnr = eff.GetXaxis().GetNbins()/float(eff.GetN()+1)
            for i, n in enumerate(nminus1s):
                eff.GetXaxis().SetBinLabel(int((i+0.5)*bnr), pretty.get(n,n))
            eff.GetXaxis().LabelsOption('v')
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1

        lg.Draw()
        ps.save(name)
        print(name, pretty_name)

# ************************************************************************
# Main
# ************************************************************************
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='submits limits for a given mass')
    parser.add_argument('--tag',help='Tagged directory name for input histograms and output plots',default='')
    parser.add_argument('--stats',help='Print stats',default=False)
    parser.add_argument('--lumi',help='Integrated luminosity in pb-1',default='1000')
    parser.add_argument('--do_tight',help='Do tight cuts',default=False)
    args = parser.parse_args()
    
    tag = args.tag
    if not os.path.exists('plots/%s'%tag):
        raise ValueError('Tagged directory %s does not exist!'%tag)

    draw_test(tag, args.stats, float(args.lumi), args.do_tight)

