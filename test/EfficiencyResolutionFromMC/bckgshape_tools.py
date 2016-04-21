#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
ps = plot_saver('plots/shapesyst')
int_lumi = 2800.

#************************************************************************************
# Helper functions
#************************************************************************************

def make_mc_hist(do_non_dy,doPI,int_lumi,low,high,rebin,hists_dir):
    '''
    - Make total MC histogram
    - inputs are [MCSamples.sample]
    - output is a TH1
    '''
    dy_samples = [dy50to120, dy120to200, dy200to400, dy400to800, dy800to1400, dy1400to2300, dy2300to3500, dy3500to4500, dy4500to6000]
    non_dy_samples = [ttbar_pow, tWtop, tWantitop, ww_incl, wz, zz_incl]
    hists_dy = []
    for sample in dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        d = f.Our2012MuonsPlusMuonsMinusHistos
        w = sample.partial_weight * int_lumi
        print sample.name, 'w = ', w
        #h = d.Get('DileptonMass').Clone('%s' % sample.name)
        h = d.Get('DimuonMassVertexConstrained').Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        #h.Draw()
        #ps.save('rawmass_%s' % sample.name)
        hists_dy.append(h)
    hists_nody = []
    for sample in non_dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        d = f.Our2012MuonsPlusMuonsMinusHistos
        w = sample.partial_weight * int_lumi
        print sample.name, 'w = ', w
        h = d.Get('DimuonMassVertexConstrained').Clone('%s' % sample.name)
        #h = d.Get('DileptonMass').Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        #h.Draw()
        #ps.save('rawmass_%s' % sample.name)
        hists_nody.append(h)
    hists = []
    # add together DY separately in case of PI background addition
    histDY = hists_dy[0].Clone('histDY')
    for j in xrange(1, len(hists_dy)):
        histDY.Add(hists_dy[j])
    histDY.Draw()
    if doPI:
        histDY = addPI(histDY,rebin,low,high)
    # now add all backgrounds
    if do_non_dy:
        hists = hists + hists_nody
        hists.append(histDY)
    else:
        hists.append(histDY)
    histMC = hists[0].Clone('histMC')
    for j in xrange(1, len(hists)):
        histMC.Add(hists[j])
    histMC.SetTitle('')
    histMC.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    histMC.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin,int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    return histMC

def addPI(histDY,rebin,low,high):
    '''
    - Multiply input DY histgram by Photon Induced background
      cross section ratio, R = (DY+PI)/DY
    - Returns TH1
    '''
    PI = ROOT.TF1('pi','[0] + [1]*x + [2]*x*x',low,high)
    PI.SetParameters(1.044e+00,2.200e-05,1.053e-08) # From Dimitri
    histPI = hist_it(PI,rebin,low,high)
    histDYPI = ROOT.TH1F('Hist', ';m(#mu^{+}#mu^{-}) [GeV];Events', 20000, 0, 20000)
    histDYPI.Rebin(rebin)
    histDYPI.GetXaxis().SetRangeUser(low, high)
    nBins = histPI.GetNbinsX()
    for i in xrange(1,nBins):
        dy = histDY.GetBinContent(i)
        dye = histDY.GetBinError(i)
        pi = histPI.GetBinContent(i)
        histDYPI.SetBinContent(i,dy*pi)
        histDYPI.SetBinError(i,dye*1.5)
    histDY.SetStats(0)
    histDY.Draw()
    histDYPI.Draw('same')
    histDYPI.SetLineColor(ROOT.kOrange+1)
    histDYPI.SetLineWidth(2)
    histDY.SetLineWidth(2)
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.AddEntry(histDY,"Nominal Background","L")
    lg.AddEntry(histDYPI,"Nominal+PI Background","L")
    lg.Draw()
    ps.save('hist_PI_overlay',log=True)
    return histDYPI

def make_pm_func(nominal, unc,uncName,fitlow,fithigh):
    '''
    - Takes in nominal background pdf and uncertainty function 
    - Returns TF1s of plus, minus background pdf shapes
    - Need to test to make naming consistent...
      At the moment the pointer name and the TF1 name need to be
      the same. Ex: func = ROOT.TF1('func',...)
    '''
    plus = ROOT.TF1('plus','(2-%s)*nominal'%uncName,fitlow,fithigh)
    minus = ROOT.TF1('minus','%s*nominal'%uncName,fitlow,fithigh)
    #plus.Draw()
    #minus.Draw('same')
    #ps.save('pm_draw_test_%s'%uncName,log=True)
    return plus, minus

def plot_fit(fit,hist,uncName,extra,low,high,fitlow,fithigh,rebin):
    '''
    - Plotting function
    - Takes in a TF1 and the histogram it was fitted to
      and plots them together. Saves png and root in
      in linear and log formats.
    - Also plots the residual of the nominal fit to the
      total MC histogram.
    '''
    # Plot Fit of Nominal vs. MC histogram
    hist.Draw()
    fit.SetLineColor(ROOT.kBlue)
    ps.c.Update()
    s = hist.GetListOfFunctions().FindObject("stats")
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(10)
    s.SetOptFit(111)
    s.Draw()
    ps.save('%s'%uncName ,log=True)
    # Plot Residual of Nominal
    res_hist(fit,hist,uncName,rebin,low,high,fitlow,fithigh)

def plot_two(f1,f1Name,hist1,f2,f2Name,hist2,uncName,extra):
    '''
    - Plotting function
      takes in 2 TF1s f1,f2 background pdf shapes
    - Returns nothing but saves a pdf, png, and root file in linear
      and log formats
    '''
    #ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetOptFit(0111)
    # Set Stats for Nominal
    #histMC.SetLineWidth(0)
    #histMC.SetLineColor(ROOT.kWhite)
    #histMC.SetTitle('Nominal Shape')
    hist1.SetName('%s'%f1Name)
    hist1.Draw()
    ps.c.Update()
    s1 = hist1.GetListOfFunctions().FindObject("stats")
    s1.SetX1NDC(0.73)
    s1.SetY1NDC(0.75)
    s1.SetOptStat(1000000001)
    s1.SetOptFit(0111)
    s1.SetName('%s'%f1Name)
    ps.c.Update()
    hist2.SetName('%s'%f2Name)
    hist2.Draw()
    ps.c.Update()
    s2 = hist2.GetListOfFunctions().FindObject("stats")
    s2.SetX1NDC(0.43)
    s2.SetX2NDC(0.70)
    s2.SetY1NDC(0.75)
    s2.SetOptStat(1000000001)
    s2.SetOptFit(0111)
    s2.SetName('%s'%f2Name)
    ps.c.Update()
    # Draw Minus
    f1.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    f1.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f1.SetLineColor(ROOT.kBlack)
    f1.SetTitle('')
    f1.Draw()
    # Draw Nominal
    f2.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f2.SetLineColor(ROOT.kOrange+1)
    f2.SetTitle('')
    s1.Draw('same')
    s2.Draw('same')
    s1.SetName('%s'%f1Name)
    s2.SetName('%s'%f2Name)
    f2.Draw("same")
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.AddEntry(f2,"%s"%f2Name,"L")
    lg.AddEntry(f1,"%s"%f1Name,"L")
    lg.Draw()
    # Save
    ps.save('%s_%s'%(uncName,extra),log=True)

def plot_three(f1,f1Name,f2,f2Name,f3,f3Name,histMC,uncName,extra=''):
    '''
    - Plotting function
      takes in 3 TF1s f1,f2,f3 and plots them
    - Returns nothing but saves a pdf, png, and root file in linear
      and log formats
    '''
    #ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetOptFit(0111)
    # Set Stats for Nominal
    #histMC.SetLineWidth(0)
    #histMC.SetLineColor(ROOT.kWhite)
    #histMC.SetTitle('Nominal Shape')
    histMC.Draw()
    ps.c.Update()
    s = histMC.GetListOfFunctions().FindObject("stats")
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(1000000001)
    s.SetOptFit(0111)
    ps.c.Update()
    # Draw f1
    f1.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    f1.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f1.SetLineColor(ROOT.kOrange+1)
    f1.SetTitle('')
    f1.Draw()
    # Draw f2
    f2.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f2.SetLineColor(ROOT.kBlue)
    f2.SetTitle('')
    f2.Draw("same")
    # Draw f3
    f3.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f3.SetLineColor(ROOT.kBlack)
    f3.SetTitle('')
    s.Draw('same')
    f3.Draw("same")
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.AddEntry(f1,f1Name,"L")
    lg.AddEntry(f2,f2Name,"L")
    lg.AddEntry(f3,f3Name,"L")
    lg.Draw()
    # Save
    name = uncName if not extra else '%s_%s'%(uncName,extra)
    ps.save('%s'%name,log=True,pdf=True)

def fit_hist_PI(hist, name,low,high,fitlow, fithigh,rebin):
    '''
    - Fit input histogram to background pdf tuned to 
      the DY+PI background shape
    - Returns a TF1
      (Redundant code from fitdymass)
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    fit.SetParNames("a", "b", "c", "d", "k")
    fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x)*x**([5])', fitlow, fithigh)
    fit.SetParNames("a", "b", "c", "d","e", "k")
    #fit.SetParameters(24, -5E-4, -5E-8, -5E-12,-5E-18, -4.5)
    # a
    fit.SetParLimits(0,-50,50)
    # b 
    fit.SetParLimits(1,-1E-1,-1E-25)
    # c
    fit.SetParLimits(2,-1E-1,-1E-25)
    # d 
    fit.SetParLimits(3,-1E-1,-1E-25)
    # e
    fit.SetParLimits(4,-1E-1,-1E-25)
    # k
    fit.SetParLimits(5,-50,50)
    '''
    hist.Fit(fit,'REMV')
    plot_fit(fit,hist,name,'',low,high,fitlow,fithigh,rebin)
    return fit

def fit_hist(hist, name,low,high,fitlow, fithigh,rebin):
    '''
    - Fit input histogram to background pdf
    - Returns a TF1
    - (Redundant code from fitdymass)
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    fit.SetParNames("a", "b", "c", "d", "k")
    fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    fit.SetParLimits(4,-10,1)
    hist.Fit(fit,'REMV')
    plot_fit(fit,hist,name,'',low,high,fitlow,fithigh,rebin)
    return fit

def hist_it(func, rebin,fitlow, fithigh):
    '''
    - Histogram-ize input function
    - Returns a TH1
    '''
    hist = ROOT.TH1F('Hist', ';m(#mu^{+}#mu^{-}) [GeV];Events', 20000, 0, 20000)
    hist.Rebin(rebin)
    hist.GetXaxis().SetRangeUser(fitlow, fithigh)
    for i in xrange(1,hist.GetNbinsX()+1):
        xlo = hist.GetXaxis().GetBinLowEdge(i)
        xhi = hist.GetXaxis().GetBinLowEdge(i+1)
        if xlo >= fitlow and xhi <= fithigh:
            integ = func.Integral(xlo,xhi) / (xhi - xlo) 
            hist.SetBinContent(i, integ)
    return hist

def fit_syst_func(nominal, histMC, unc, uncName, rebin,low,high,fitlow, fithigh,extra):
    '''
    - Take plus/minus function, create a histogram identical to input MC histogram
      out of it and fit it with the background pdf
    - uncName needs to be consistent with bckgshape_syst module (e.g. unc50)
    - Returns 2 TF1s (doesn't do this)
    '''
    plus, minus = make_pm_func(nominal,unc,uncName,fitlow,fithigh)
    #plot_fit(nominal,histMC,uncName,'func',low,high,fitlow,fithigh,rebin)
    plot_three(minus,extra+'% Decrease',plus,extra+'% Increase',nominal,'Nominal',histMC,uncName,'func')
    res_hist(plus,histMC,uncName+'plus',rebin,low,high,fitlow,fithigh,resmin=-1,resmax=1)
    res_hist(minus,histMC,uncName+'minus',rebin,low,high,fitlow,fithigh,resmin=-1,resmax=1)
    #return plus,minus

def res_hist(fit,hist,uncName,rebin,low,high,fitlow,fithigh,resmin=-0.25,resmax=0.25):
    '''
    - Create Residual histogram (fit-hist)/hist
    - Inputs are TF1 fit to a TH1 and the TH1
    - Returns nothing but saves the residual histogram of the fit and 
      the input data to the fit
    '''
    hres = ROOT.TH1F('hres', ';m(#mu^{+}#mu^{-}) [GeV];(fit-hist)/hist', 20000, 0, 20000)
    hres.Rebin(rebin)
    xax = hres.GetXaxis()
    xax.SetRangeUser(low, high)
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= low and xhi <= high:
            res = fit.Integral(xlo, xhi)/(xhi-xlo) - hist.GetBinContent(i)
            if hist.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/hist.GetBinContent(i))
                hres.SetBinError(i, hist.GetBinError(i)/hist.GetBinContent(i))
    hres.SetMinimum(resmin)
    hres.SetMaximum(resmax)
    hres.Draw('e')
    l1 = ROOT.TLine(fitlow, 0., fithigh,  0.)
    l1.Draw()
    ps.save('res_%s'%uncName, log=False)

def make_pm_hist(histMC, uncHist,name,rebin,low,high):
    '''
    - Takes in two histograms 1 nominal and 2 shape systematic
    - Returns two TH1s plus, minus shapes
    '''
    plusHist = ROOT.TH1F('hres_%s' %name, ';m(#mu^{+}#mu^{-}) [GeV];(fit-hist)/hist', 20000, 0, 20000)
    plusHist.Rebin(rebin)
    plusHist.GetXaxis().SetRangeUser(low, high)
    minusHist = plusHist.Clone('minusHist')
    # loop on bins and multiply them together
    for i in xrange(1,histMC.GetNbinsX()+1):
        plus = (2-uncHist.GetBinContent(i))*histMC.GetBinContent(i)
        minus = uncHist.GetBinContent(i)*histMC.GetBinContent(i)
        plusHist.SetBinContent(i, plus)
        minusHist.SetBinContent(i, minus)
    return plusHist, minusHist
        
def fit_syst_hist(histMC,nominal,unc,uncName,rebin,fitlow,fithigh,low,high):
    '''
    - Take input MC histogram, apply bckgshape_syst +/- functions to it then fit it
      with the background pdf. Returns two TF1s.
    - Alternative method to obtain +/- fit functions.
      Not used.
    '''
    # Histogramize bckgshape_syst.unc
    uncHist = hist_it(unc,rebin,fitlow,fithigh)
    # Multiply hist_syst 
    histPlus,histMinus = make_pm_hist(histMC,uncHist,uncName,rebin,low,high)
    plusFit = fit_hist(histPlus,'plusFit_%s'%uncName,low,high,fitlow,fithigh,rebin)
    minusFit = fit_hist(histPlus,'minusFit_%s'%uncName,low,high,fitlow,fithigh,rebin)
    plot_three(plusFit,'Plus',minusFit,'Minus',nominal,'Nominal',uncName,'hist')

#************************************************************************************
# Running the code interactively to test new modules
#************************************************************************************

if __name__=='__main__':
    '''
    - Running module as stand alone code
    - May not work as expected if tested from
      fitMCMassSpectra.py
    '''
    import bckgshape_syst
    #hists_dir = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/DataMCSpectraComparision/mc/76X_v2/' # should be here oops
    hists_dir = '/afs/cern.ch/work/c/cschnaib/DataMCSpectraComparision/mc/76X_v2/'

    set_zp2mu_style()
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.TH1.AddDirectory(0)
    ps = plot_saver('plots/testTools')

    rebin = 40
    low = 40
    fitlow = 400
    fithigh = 5500
    high = 5500
    # Run2015(B+C+D)
    int_lumi = 2800.

    # Make nominal MC histogram
    histMC = make_mc_hist(True,False,int_lumi,low,high,rebin,hists_dir)
    # Make MC histogram with PI background
    histMC = make_mc_hist(True,True,int_lumi,low,high,rebin,hists_dir)
    # Make nominal MC fit
    nominal = fit_hist(histMC,'nominal',rebin,fitlow,fithigh,rebin)
    #fit_syst_func(nominal,histMC,bckgshape_syst.unc1,'unc1',fitlow,fithigh,low,high)
    #fit_syst_func(nominal,histMC,bckgshape_syst.unc2,'unc2',fitlow,fithigh,low,high)
    #fit_syst_func(nominal,histMC,bckgshape_syst.unc3,'unc3',fitlow,fithigh,low,high)
    fit_syst_hist(histMC,nominal,bckgshape_syst.unc1,'unc1',rebin,fitlow,fithigh,low,high)
