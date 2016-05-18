#!/usr/bin/env python

# (py draw.py >! plots/nminus1effs/out.draw) && tlp plots/nminus1effs

# Import nm1entry.py
# - from MCSamples and roottools import *
# - nm1entry class
# - nminus1s list
# - pretty dictionary
# - styles dictionary
from nm1entry import *

def draw_mass_test(tag, printStats, lumi, do_tight=False):

    # Specific stylings for this script on top of zp2mu_style
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ROOT.gStyle.SetTitleX(0.12)
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

    refXS = dy50to120.cross_section
    refN = dy50to120.nevents

    samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]
    # All MC samples
    # lumi
    mc_samples = [nm1entry(sample,False,lumi) for sample in samples]
    #for mc_sample in mc_samples:
    #    exec '%s = mc_sample' % mc_sample.name

    #bin_width = 20
    #maxX = 2500
    #minX = 60
    #nBins = (maxX-minX)/bin_width
    #mass_range = []
    #for i in range(3,nBins):
    #    ibin = i*bin_width
    #    mass_range.append(ibin)
    #print mass_range
    mass_range = [60,120,180,240,320,500,1000,2500]

    to_use = {
    #   'sample':[MC,Data],
        'NoPt':[mc_samples,data],
        'NoDB':[mc_samples,data],
        'NoIso':[mc_samples,data],
        'NoTkLayers':[mc_samples,data],
        'NoPxHits':[mc_samples,data],
        'NoMuHits':[mc_samples,data],
        'NoMuMatch':[mc_samples,data],
        'NoVtxProb':[mc_samples,data],
        'NoB2B':[mc_samples,data],
        'NoDptPt':[mc_samples,data],
        'NoTrgMtch':[mc_samples,data],
        }

    yrange = {
    #   'sample':    (ymin,ymax),
        'NoPt':      (0.00,1.01),
        'NoDB':      (0.95,1.001),
        'NoIso':     (0.60,1.01),
        'NoTkLayers':(0.95,1.001),
        'NoPxHits':  (0.80,1.001),
        'NoMuHits':  (0.80,1.001),
        'NoMuMatch': (0.80,1.005),
        'NoVtxProb': (0.90,1.001),
        'NoB2B':     (0.80,1.001),
        'NoDptPt':   (0.95,1.001),
        'NoTrgMtch': (0.90,1.001),
        }
    #global_ymin = 0.
    global_ymin = None


    ROOT.gStyle.SetTitleX(0.25)
    ROOT.gStyle.SetTitleY(0.50)

    for nminus1 in nminus1s:
        pretty_name = pretty[nminus1]
        lg = ROOT.TLegend(0.25, 0.21, 0.91, 0.44)
        lg.SetTextSize(0.03)
        lg.SetFillColor(0)
        lg.SetBorderSize(1)
        same = 'A'
        effs = []
        for entry in to_use[nminus1]: #,mass_range 
            l = len(mass_range)-1
            nminus1_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            nminus1_den = ROOT.TH1F('den', '', l, array('f',mass_range))
            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                if printStats: table(entry,nminus1,mass_range)
                color, fill = styles[entry.name]
                hnum = entry.histos['NoNo']
                hden = entry.histos[nminus1]
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                    nminus1_num.SetBinContent(mbin+1, num)
                    nminus1_den.SetBinContent(mbin+1, den)
                eff,p,epl,eph = binomial_divide(nminus1_num, nminus1_den)
            else:
                if printStats: table_wald(mc,nminus1,mass_range)
                p_hats = []
                errsW = []
                x = []
                ex = []
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    numTot = 0
                    denTot = 0
                    err2sum = 0
                    numsW = []
                    densW = []
                    err2s = []
                    for i,mc in enumerate(entry):
                        hnum = mc.histos['NoNo']
                        hden = mc.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        numInt = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                        denInt = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                        if numInt<=denInt and denInt>0:
                            p_hat_mc = float(numInt)/denInt
                            err2 = p_hat_mc*(1-p_hat_mc)/denInt
                        elif numInt>denInt:
                            print "numInt>denInt"
                            print "NM1, sample, mlow, mhigh, num, den"
                            print nminus1, mc.name, mlow, mhigh, numInt, denInt
                            print
                            p_hat_mc = 1
                            err2 = 0
                        else:
                            if numInt!=0 or denInt!=0:
                                print "ELSE"
                                print "NM1, sample, mlow, mhigh, num, den"
                                print nminus1, mc.name, mlow, mhigh, numInt, denInt
                                print
                            p_hat_mc = 0
                            err2 = 0
                        numTot = numTot + numInt*mc.partial_weight
                        denTot = denTot + denInt*mc.partial_weight
                        numsW.append(numInt*mc.partial_weight)
                        densW.append(denInt*mc.partial_weight)
                        err2s.append(err2)
                    p_hat = float(numTot)/denTot
                    err2sum = sum( ((m/denTot)**2 * e2) for m,e2 in zip(densW,err2s))
                    if err2sum<0:
                        err2sum = 0
                    elif p_hat==0:
                        err2sum = 0
                    p_hats.append(p_hat)
                    errsW.append(err2sum**0.5)
                    x.append(nminus1_num.GetXaxis().GetBinCenter(mbin+1))
                    ex.append(nminus1_num.GetXaxis().GetBinWidth(mbin+1)/2)
                eff = ROOT.TGraphAsymmErrors(len(x), *[array('d',obj) for obj in (x,p_hats,ex,ex,errsW,errsW)])
            eff.SetTitle(pretty_name)
            ymin, ymax = yrange[nminus1]
            eff.GetYaxis().SetRangeUser(global_ymin if global_ymin is not None else ymin, ymax)
            eff.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
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
                if len(entry)==len(mc_samples):
                    draw = '2'
                    eff.SetLineColor(ROOT.kGreen+2)
                    eff.SetFillColor(ROOT.kGreen+2)
                    eff.SetFillStyle(1001)
                    lg.AddEntry(eff,'All Simulation','LF')
                else:
                    raise ValueError("Set line+fill color, fill Style, and set legend entry for MC entry list")
            #lg.AddEntry(eff, pretty.get(entry.name, entry.name), 'LF')
            draw += same
            eff.Draw(draw)
            effs.append(eff)
            same = ' same'
            outfile.cd()
            eff.Write("arp%d"%iarp)
            iarp+=1
        # end for entry in to_use[name]: # entry is a specific sample
        lg.Draw()
        name = nminus1+'_mass'
        ps.save(name)
        print(nminus1, pretty_name, name)
    # end for name, mass_range in mass_bins:

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

    # cjsbad - add a cmd-line option for dy only, qcd, zprime, etc
    tag = args.tag
    if not os.path.exists('plots/%s'%tag):
        raise ValueError('Tagged directory %s does not exist!'%tag)

    draw_mass_test(tag, args.stats, float(args.lumi), args.do_tight)

