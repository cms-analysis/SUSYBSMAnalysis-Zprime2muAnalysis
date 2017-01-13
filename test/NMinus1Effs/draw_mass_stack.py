#!/usr/bin/env python

# Import nm1entry.py
# - from MCSamples and roottools import *
# - nm1entry class
# - nminus1s list
# - pretty dictionary
# - styles dictionary
from nm1entry import *
ROOT.gROOT.SetBatch(True)

def draw_mass_stack_test(tag, printStats, lumi, do_tight=False):

    # Specific stylings for this script on top of zp2mu_style
    #ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.02)
    #ROOT.gStyle.SetTitleX(0.12)
    #ROOT.gStyle.SetTitleH(0.07)
    #ROOT.gStyle.SetTitleX(0.25)
    #ROOT.gStyle.SetTitleY(0.50)
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
    ps.c.SetTopMargin(0.075)
    ps.c.SetBottomMargin(0.1)
    ps.c.SetRightMargin(0.05)
    #ps.c.SetLeftMargin(0.03)
    ps.c.SetLogy(1)

    # cjsbad
    # fix the data entry
    data = nm1entry('data', True, lumi)

    refXS = dy50to120.cross_section
    refN = dy50to120.nevents
    #print lumi, refN/refXS

    # All MC samples
    samples = [dy50to120,dy120to200,dy200to400,dy400to800,dy800to1400,dy1400to2300,dy2300to3500,dy3500to4500,dy4500to6000,ttbar_pow,ww_incl,zz_incl,wz,tWtop,tWantitop,qcd80to120,qcd120to170,qcd170to300,qcd300to470,qcd470to600,qcd600to800,qcd800to1000,qcd1000to1400,qcd1400to1800,qcd1800to2400,qcd2400to3200,qcd3200,wjets]
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


    print nminus1s
    for nminus1 in nminus1s:
        pretty_name = pretty[nminus1]
        print nminus1, pretty_name
        lg = ROOT.TLegend(0.45, 0.55, 0.91, 0.9)
        lg.SetTextSize(0.03)
        lg.SetFillColor(0)
        lg.SetBorderSize(1)
        
        same = 'A'
        effs = []

        stack_num = ROOT.THStack('hs','')
        stack_den = ROOT.THStack('hs','')



        for entry in to_use[nminus1]: #,mass_range 

            l = len(mass_range)-1
            data_num = ROOT.TH1F('num', '', l, array('f',mass_range))
            data_den = ROOT.TH1F('den', '', l, array('f',mass_range))

            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                if printStats: table(entry,nminus1, mass_range)
                color, fill = styles[entry.name]
                hnum = entry.histos['NoNo']
                hden = entry.histos[nminus1]
                for mbin in range(len(mass_range)):
                    if mbin == (len(mass_range)-1): continue
                    mlow = mass_range[mbin]
                    mhigh = mass_range[mbin+1]
                    num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                    den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                    data_num.SetBinContent(mbin+1, num)
                    data_den.SetBinContent(mbin+1, den)
            else:
                if printStats: table_wald(mc,nminus1,mass_range)
                for i,mc in enumerate(entry):
                    nminus1_den_MC = ROOT.TH1F('den', '', l, array('f',mass_range))
                    nminus1_num_MC = ROOT.TH1F('num', '', l, array('f',mass_range))
                    for mbin in range(len(mass_range)):
                        if mbin == (len(mass_range)-1): continue
                        hnum = mc.histos['NoNo']
                        hden = mc.histos[nminus1]
                        mlow = mass_range[mbin]
                        mhigh = mass_range[mbin+1]
                        num = get_integral(hnum, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                        den = get_integral(hden, mlow, mhigh, integral_only=True, include_last_bin=False, nm1=True)
                        nminus1_num_MC.SetBinContent(mbin+1, num)
                        nminus1_den_MC.SetBinContent(mbin+1, den)
                    color, fill = styles[mc.name]
                    nminus1_num_MC.SetLineColor(color)
                    nminus1_num_MC.SetFillColor(color)
                    nminus1_num_MC.SetMinimum(1E-7)
                    #nminus1_num_MC.Scale(mc.partial_weight*refN/refXS)
                    nminus1_num_MC.Scale(mc.partial_weight*lumi)
                    stack_num.Add(nminus1_num_MC)
                    nminus1_den_MC.SetLineColor(color)
                    nminus1_den_MC.SetFillColor(color)
                    #nminus1_den_MC.Scale(mc.partial_weight*refN/refXS)
                    nminus1_den_MC.Scale(mc.partial_weight*lumi)
                    stack_den.Add(nminus1_den_MC)
                    #lg.AddEntry(nminus1_den_MC,pretty.get(mc.name,mc.name),"F")
                    if mc.name == 'dy50to120':
                        lg.AddEntry(nminus1_den_MC,"Drell-Yan",'F' )
                    elif mc.name == 'ttbar_pow':
                        lg.AddEntry(nminus1_den_MC,"t#bar{t}","F")
                    elif mc.name == 'ww_incl':
                        lg.AddEntry(nminus1_den_MC,"DiBoson","F")
                    elif mc.name == 'qcd80to120':
                        lg.AddEntry(nminus1_den_MC,"QCD & W+jets","F")
                    elif mc.name == 'tWtop':
                        lg.AddEntry(nminus1_den_MC,"Single Top","F")
            if not isinstance(entry,(list,tuple)) and 'data' in entry.name:
                draw = 'P'
                data_num.SetLineColor(color)
                data_num.SetMarkerStyle(20)
                data_num.SetMarkerSize(1.05)
                data_num.SetMarkerColor(color)
                data_den.SetLineColor(color)
                data_den.SetMarkerStyle(20)
                data_den.SetMarkerSize(1.05)
                data_den.SetMarkerColor(color)
                #lg.AddEntry(data, pretty.get(entry.name, entry.name) % (lumi/1000.), 'LP')
                lg.AddEntry(data_den, pretty.get(entry.name, entry.name) % (entry.lumi/1000.), 'LP')
        #stack_num.SetMinimum(1E-7)
        stack_num.Draw("hist")
        data_num.Draw("pe1same")
        lg.Draw()
        outfile.cd()
        #stack_num.SetMinimum(0.1)
        stack_num.SetTitle("All Selection Applied")
        stack_num.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
        #stack_num.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        stack_num.GetYaxis().SetTitle("Events")
        numName = nminus1+'_stack_mass_num'
        ps.save(numName)
        stack_den.Draw("hist")
        data_den.Draw("pe1same")
        lg.Draw()
        outfile.cd()
        stack_den.SetMinimum(0.1)
        stack_den.SetTitle("N - ("+pretty_name+")")
        stack_den.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
        #stack_den.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        stack_den.GetYaxis().SetTitle("Events")
        denName = nminus1+'_stack_mass_den'
        ps.save(denName)
        print(nminus1, pretty_name, denName, numName)
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

    tag = args.tag
    if not os.path.exists('plots/%s'%tag):
        raise ValueError('Tagged directory %s does not exist!'%tag)

    draw_mass_stack_test(tag, args.stats, float(args.lumi), args.do_tight)

