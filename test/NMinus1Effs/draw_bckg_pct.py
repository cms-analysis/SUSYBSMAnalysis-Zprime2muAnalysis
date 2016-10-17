#!/usr/bin/env python

# Import nm1entry.py
# - from MCSamples and roottools import *
# - nm1entry class
# - nminus1s list
# - pretty dictionary
# - styles dictionary
from nm1entry import *
ROOT.gROOT.SetBatch(True)

def draw_bckg_pct_test(tag, printStats, lumi):

    # Specific stylings for this script on top of zp2mu_style
    ROOT.gStyle.SetPadRightMargin(0.03)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.07)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.TH1.AddDirectory(0)
    # output (find a cmd-line way for psn)
    outfile = ROOT.TFile("test.root","recreate")

    # 'plots' = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/NMinus1Effs/plots/TAG/tag'
    # TAG = MC/Data production TAG
    # tag = sub tag for each version of plots made with a single production TAG
    psn = 'plots'
    if tag:
        psn = 'plots/%s'%tag
    ps = plot_saver(psn, size=(600,600), log=True, pdf=True, name='bckg_plot')

    samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]
    mc_samples = [nm1entry(sample,False,lumi) for sample in reversed(samples)]
    for mc_sample in mc_samples:
        exec '%s = mc_sample' % mc_sample.name

    # Include 'NoNo' in nminus1s list
    nminus1s.append('NoNo')

    to_use = {
        'NoPt':[mc_samples],
        'NoDB':[mc_samples],
        'NoIso':[mc_samples],
        'NoTkLayers':[mc_samples],
        'NoPxHits':[mc_samples],
        'NoMuHits':[mc_samples],
        'NoMuMatch':[mc_samples],
        'NoVtxProb':[mc_samples],
        'NoB2B':[mc_samples],
        'NoDptPt':[mc_samples],
        'NoTrgMtch':[mc_samples],
        'NoNo':[mc_samples],
    }

    mass_range = [60,80,100,120]
    maxX = 5000
    bin_width = 50
    nBins = maxX/bin_width
    for i in xrange(3,nBins):
        ibin = i*bin_width
        mass_range.append(ibin)
    nmc = len(mc_samples)
    nnm1 = len(nminus1s)
    nbins = len(mass_range)-1
    sdata = (nnm1,nbins,nmc)
    stotals = (nnm1,nbins)
    data = np.zeros(sdata)
    totals = np.zeros(stotals)

    #
    # There may be a way to do this in one pass through the input histos, stick to two for now
    # 

    # - Make data
    #   - data matrix 
    #     3 dimensions (nminus1,mass_bin_lowEdge,mc_sample)
    #   - totals matrix
    #     2 dimensions (nminus1,mass_bin_lowEdge)
    for nm1,nminus1 in enumerate(nminus1s):
        for entry in to_use[nminus1]: # single item (mc_samples)
            for ibin in range(len(mass_range)):
                if ibin == nbins: continue
                mlow = mass_range[ibin]
                mhigh = mass_range[ibin+1]
                total = 0
                for imc,mc in enumerate(entry):
                    hmc = mc.histos[nminus1]
                    count = get_integral(hmc,mlow,mhigh,integral_only=True,include_last_bin=False, nm1=True)
                    count = count*sample.partial_weight
                    data.itemset((nm1,ibin,imc),count)
                    total = total+count
                totals.itemset((nm1,ibin),total)

    # Make the Histograms
    for nm1,nminus1 in enumerate(nminus1s):
        pretty_name = pretty[nminus1]
        lg = ROOT.TLegend(0.65,0.65,0.95,0.85)
        stack = ROOT.THStack('hs','')
        for entry in to_use[nminus1]:
            for mc,sample in enumerate(entry):
                mc_hist = ROOT.TH1F('mcHist','',nbins,array('f',mass_range))
                mc_hist.SetMaximum(1.)
                for ibin in range(len(mass_range)):
                    if ibin == nbins: continue
                    binContent = data[nm1,ibin,mc]/totals[nm1,ibin]
                    mc_hist.SetBinContent(ibin,binContent)
                color,fill = styles[sample.name]
                mc_hist.SetLineColor(color)
                mc_hist.SetFillColor(color)
                stack.Add(mc_hist)
                if sample.name == 'dy50to120':
                    lg.AddEntry(mc_hist,"Drell-Yan",'F' )
                elif sample.name == 'ttbar_pow':
                    lg.AddEntry(mc_hist,"t#bar{t}","F")
                elif sample.name == 'ww_incl':
                    lg.AddEntry(mc_hist,"DiBoson","F")
                elif sample.name == 'qcd80to120':
                    lg.AddEntry(mc_hist,"QCD & W+jets","F")
                elif sample.name == 'tWtop':
                    lg.AddEntry(mc_hist,"Single Top","F")
        stack.SetMaximum(1.)
        stack.SetMinimum(1.E-4)
        stack.Draw('hist')
        stack.SetTitle(pretty_name+' Background Percentage')
        stack.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
        stack.GetYaxis().SetTitle('percentage/%i GeV'%bin_width)
        stack.GetYaxis().SetTitleOffset(1.4)
        lg.Draw()
        name = nminus1+'_bckg_pct'
        print(nminus1,pretty_name,name)
        ps.save(name,pdf=True,log=True)

# ************************************************************************
# Main
# ************************************************************************
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='submits limits for a given mass')
    parser.add_argument('--tag',help='Tagged directory name for input histograms and output plots',default='')
    parser.add_argument('--stats',help='Print stats',default=False)
    parser.add_argument('--lumi',help='Integrated luminosity in pb-1',default='1000')
    args = parser.parse_args()

    tag = args.tag
    if not os.path.exists('plots/%s'%tag):
        raise ValueError('Tagged directory %s does not exist!'%tag)

    draw_bckg_pct_test(tag, args.stats, float(args.lumi))

