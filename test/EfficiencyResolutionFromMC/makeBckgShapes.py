#!/usr/bin/env python

##############################
# Histogram making functions #
##############################

def make_mc_hist(do_non_dy,doPI,int_lumi,low,high,rebin,hists_dir,rootDir,rootHist,title,cat):
    '''
    - Make total MC histogram
    - inputs are [MCSamples.sample]
    - output is a TH1
    '''
    dy_samples = [dy50to120, dy120to200, dy200to400, dy400to800, dy800to1400, dy1400to2300, dy2300to3500, dy3500to4500, dy4500to6000,dyInclusive50]
    non_dy_samples = [ ttbar_lep, tW, Wantitop, WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, WZ, ZZ]#ttbar_pow,
    #non_dy_samples = [ttbar_pow, tW, Wantitop, WWinclusive, WZ, ZZ]
    hists_dy = []
    for sample in dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        w = sample.partial_weight * int_lumi
        #print sample.name, 'w = ', w
        #h = d.Get('DileptonMass').Clone('%s' % sample.name)
        h = f.Get(rootDir).Get(rootHist).Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        hists_dy.append(h)
    hists_nody = []
    for sample in non_dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        w = sample.partial_weight * int_lumi
        #print sample.name, 'w = ', w
        h = f.Get(rootDir).Get(rootHist).Clone('%s' % sample.name)
        #h = d.Get('DileptonMass').Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        hists_nody.append(h)
    hists = []
    # add together DY separately in case of PI background addition
    histDY = hists_dy[0].Clone('histDY')
    for j in xrange(1, len(hists_dy)):
        histDY.Add(hists_dy[j])
    histDY.Draw()
    if doPI:
        histDY = addPI(histDY,rebin,low,high,int_lumi,title,cat)
    # now add all backgrounds
    if do_non_dy:
        hists = hists + hists_nody
        hists.append(histDY)
    else:
        hists.append(histDY)
    histMC = hists[0].Clone('histMC')
    for j in xrange(1, len(hists)):
        histMC.Add(hists[j])
    histMC.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    histMC.GetXaxis().SetLabelFont(42)
    # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    histMC.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin,int_lumi/1000)) 
    histMC.GetYaxis().SetLabelFont(42)
    return histMC

def addPI(histDY,rebin,low,high,int_lumi,title,cat):
    '''
    - Multiply input DY histgram by Photon Induced background
      cross section ratio, R = (DY+PI)/DY
    - Returns TH1
    '''
    PI = ROOT.TF1('pi','[0] + [1]*x + [2]*x*x',low,high)
    #PI.SetParameters(1.044e+00,2.200e-05,1.053e-08) # From Dimitri Method 0 - Paper 2016
    PI.SetParameters(1.025e+00,7.188e-06,1.457e-09) # From Dimitri Method 1 - ICHEP 2016
    histPI = hist_it(PI,rebin,low,high)
    histDYPI = ROOT.TH1F('Hist', ';m(#mu^{+}#mu^{-}) [GeV];Events', 20000, 0, 20000)
    histDYPI.Rebin(rebin)
    histDYPI.GetXaxis().SetRangeUser(low, high)
    histDYPI.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    histDYPI.GetXaxis().SetLabelFont(42)
    # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    histDYPI.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin,int_lumi/1000)) 
    histDYPI.GetYaxis().SetLabelFont(42)
    nBins = histPI.GetNbinsX()
    for i in xrange(1,nBins):
        dy = histDY.GetBinContent(i)
        dye = histDY.GetBinError(i)
        pi = histPI.GetBinContent(i)
        histDYPI.SetBinContent(i,dy*pi)
        # why? cjsbad
        histDYPI.SetBinError(i,dye*1.5)
    histDY.SetStats(0)
    histDYPI.SetStats(0)
    histDYPI.Draw()
    histDY.Draw('same')
    histDYPI.SetLineColor(ROOT.kOrange+1)
    histDYPI.SetLineWidth(2)
    histDY.SetLineWidth(2)
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(histDY,'Nominal Background','LE')
    lg.AddEntry(histDYPI,'Nominal+PI Background','LE')
    lg.Draw()
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    ps.save('%s_PI_hist_overlay'%(cat),pdf=True,log=True)
    return histDYPI

def res_hist(fit,hist,rebin,low,high,resmin=-0.25,resmax=0.25):
    '''
    - Inputs are TF1 fit to a TH1 and the TH1
    - Returns TH1 Residual histogram (fit-hist)/hist
    '''
    hres = ROOT.TH1F('hres', ';m(#mu^{+}#mu^{-}) [GeV];(fit-hist)/hist', 20000, 0, 20000)
    hres.GetXaxis().SetLabelFont(42)
    hres.GetYaxis().SetLabelFont(42)
    hres.Rebin(rebin)
    xax = hres.GetXaxis()
    xax.SetRangeUser(low, high)
    for h in [hres]:
        h.SetMarkerStyle(2)
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
    return hres

#####################
# Fitting functions #
#####################

def fit_data(hist, name,low,high,fitlow, fithigh,rebin,title):
    '''
    - Fit input histogram to background pdf
    - Returns a TF1
    '''
    #fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    #fit.SetParNames('a', 'b', 'c', 'd', 'k')
    #fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x)*x**([3])', fitlow, fithigh)
    fit.SetParNames('a', 'b', 'c', 'k')
    fit.SetParameters(24, -5E-4, -5E-8, -4.5)
    #fit.SetParLimits(2,-1,0)
    #fit.SetParLimits(3,-1,0)
    print '\nFit Data Shape',name,'\n'
    hist.Fit(fit,'REM')
    plot_fit(fit,hist,name,low,high,fitlow,fithigh,rebin,title,data=True)
    return fit

def fit_bckg(hist, name,low,high,fitlow, fithigh,rebin,title):
    '''
    - Fit input histogram to background pdf
    - Returns a TF1
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    fit.SetParNames('a', 'b', 'c', 'd', 'k')
    fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    fit.SetParLimits(4,-10,1)
    print '\nFit Background Shape',name,'\n'
    hist.Fit(fit,'REM')
    plot_fit(fit,hist,name,low,high,fitlow,fithigh,rebin,title)
    return fit

def fit_bckg_PI(hist, name,low,high,fitlow, fithigh,rebin,title):
    '''
    - Fit input histogram to background pdf tuned to 
      the DY+PI background shape
    - Returns a TF1
    '''
    fit = ROOT.TF1('%s'%name, 'exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])', fitlow, fithigh)
    fit.SetParNames('a', 'b', 'c', 'd', 'k')
    fit.SetParameters(24, -5E-4, -5E-8, -5E-12, -4.5)
    print '\nFit Background Shape with PI',name,'\n'
    hist.Fit(fit,'REM')
    plot_fit(fit,hist,name,low,high,fitlow,fithigh,rebin,title)
    return fit

def fit_line(hist,fitName,fitlow,fithigh):
    '''
    - Fits residual histograms to a straight line
    - Returns a TF1
    '''
    fit = ROOT.TF1('%s'%fitName, 'pol0', fitlow, fithigh)
    fit.SetParNames('r')
    print '\nFit line to residual plot',fitName,'\n'
    hist.Fit(fit,'REM')
    fit.SetLineWidth(2)
    fit.SetLineColor(ROOT.kGreen+2)
    return fit

def make_data_hist(data_dir,rootDir,rootHist,rebin,low,high,title,cat,int_lumi,varBins=False):
    '''
    - Fit bckg shape to data
    '''
    fn ='ana_datamc_data.root'
    fn = data_dir + fn
    f = ROOT.TFile(fn)
    dataHist = f.Get(rootDir).Get(rootHist)
    if varBins == True:
        rebinning = array.array('d',[])
        for i in range(low,400,10):
            rebinning.append(i)
        for i in range(400,800,15):
            rebinning.append(i)
        for i in range(800,1000,20):
            rebinning.append(i)
        for i in range(1000,1200,25):
            rebinning.append(i)
        for i in range(1200,1500,50):
            rebinning.append(i)
        for i in range(1500,high+1,100):
            rebinning.append(i)
        #for i in range(1500,high+1,int(high-1500)):
        #    rebinning.append(i)
        #print rebinning
        dataHist = dataHist.Rebin(len(rebinning)-1,"%s_copy"%(dataHist.GetName()),rebinning)
    else: dataHist.Rebin(rebin)
    dataHist.GetXaxis().SetRangeUser(low, high)
    dataHist.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    dataHist.GetXaxis().SetLabelFont(42)
    # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    dataHist.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin,int_lumi/1000)) 
    #dataHist.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) 
    dataHist.GetYaxis().SetLabelFont(42)
    dataHist.Draw('hist')
    ps.save('test_dataHist_%s'%cat,log=True)
    return dataHist
    

######################
# Plotting functions #
######################

def plot_fit(fit,hist,fitName,low,high,fitlow,fithigh,rebin,title,data=False):
    '''
    - Plotting function
    - Takes in a TF1 and the histogram it was fitted to
      and plots them together. Saves png and root in
      in linear and log formats.
    - Also plots the residual of the nominal fit to the
      total MC histogram.
    '''
    # Plot Fit of Nominal vs. MC histogram
    hist.SetTitle('')
    hist.Draw()
    fit.SetLineColor(ROOT.kBlue)
    ps.c.Update()
    s = hist.GetListOfFunctions().FindObject('stats')
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(10)
    s.SetOptFit(1111)
    s.Draw()
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    if data==True: hist.SetMinimum(1E-2)
    ps.save('%s'%(fitName) ,pdf=True,log=True)
    # Plot Residual of Nominal
    hres = res_hist(fit,hist,rebin,low,high)
    line = fit_line(hres,fitName,fitlow,fithigh)
    plot_res(hres,line,fitName,title,fitlow,fithigh)

def plot_two(f1,f1Name,hist1,f2,f2Name,hist2,uncName,extra,title,cat,int_lumi):
    '''
    - Plotting function
      takes in 2 TF1s f1,f2 background pdf shapes
    - Returns nothing but saves a pdf, png, and root file in linear
      and log formats
    '''
    ROOT.gStyle.SetOptStat("n");
    # Set Stats for Nominal
    hist1.SetName('%s'%f1Name)
    hist1.Draw()
    ps.c.Update()
    s1 = hist1.GetListOfFunctions().FindObject('stats')
    s1.SetX1NDC(0.73)
    s1.SetY1NDC(0.75)
    s1.SetOptStat(0)
    s1.SetOptFit(1111)
    s1.SetName('%s'%f1Name)
    ps.c.Update()
    hist2.SetName('%s'%f2Name)
    hist2.Draw()
    ps.c.Update()
    s2 = hist2.GetListOfFunctions().FindObject('stats')
    s2.SetX1NDC(0.73)
    s2.SetY1NDC(0.48)
    s2.SetY2NDC(0.73)
    s2.SetOptStat(0)
    s2.SetOptFit(1111)
    s2.SetName('%s'%f2Name)
    ps.c.Update()
    # Draw Minus
    f1.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    f1.GetYaxis().SetLabelFont(42)
    f1.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f1.GetXaxis().SetLabelFont(42)
    f1.SetLineColor(ROOT.kBlack)
    f1.SetTitle('')
    f1.Draw()
    # Draw Nominal
    f2.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    f2.SetLineColor(ROOT.kOrange+1)
    f2.SetTitle('')
    s1.Draw('same')
    s2.Draw('same')
    f2.Draw('same')
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(f2,'%s'%f2Name,'L')
    lg.AddEntry(f1,'%s'%f1Name,'L')
    lg.Draw()
    # Draw label
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Save
    ps.save('%s_%s_%s'%(cat,uncName,extra),pdf=True,log=True)

def plot_three(hist,fitNom,fitPlus,fitMinus,uncName,cat,title,int_lumi):
    '''
    - Plotting function
      takes in 3 TF1s fitNom,fitPlus,fitMinus and plots them
    - Returns nothing but saves a pdf, png, and root file in linear
      and log formats
    '''
    # Settigns fitNom
    fitNom.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    fitNom.GetYaxis().SetLabelFont(42)
    fitNom.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    fitNom.GetXaxis().SetLabelFont(42)
    fitNom.SetLineColor(ROOT.kBlack)
    fitNom.SetTitle('')
    fitNom.GetXaxis().SetRangeUser(200,5500)
    # Settings fitPlus
    fitPlus.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    fitPlus.GetYaxis().SetLabelFont(42)
    fitPlus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    fitPlus.GetXaxis().SetLabelFont(42)
    fitPlus.SetLineColor(ROOT.kBlue)
    fitPlus.SetTitle('')
    fitPlus.GetXaxis().SetRangeUser(200,5500)
    # Settings fitMinus
    fitMinus.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) # int_lumi is in pb-1.
    fitMinus.GetYaxis().SetLabelFont(42)
    fitMinus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
    fitMinus.GetXaxis().SetLabelFont(42)
    fitMinus.SetLineColor(ROOT.kOrange+1)
    fitMinus.SetTitle('')
    fitMinus.GetXaxis().SetRangeUser(200,5500)
    # Draw Fits
    fitMinus.Draw('L')
    fitNom.Draw('Lsame')
    fitPlus.Draw('Lsame')
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    uncLabel = uncName+'%'
    lg.AddEntry(fitNom,'Nominal','L')
    lg.AddEntry(fitPlus,'+ %s'%uncLabel,'L')
    lg.AddEntry(fitMinus,'- %s'%uncLabel,'L')
    lg.Draw()
    # Draw label
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Save
    ps.save('%s_Fit_%spm'%(cat,uncName),pdf=True,log=True)

def plot_res(hres,line,fitName,title,fitlow,fithigh):
    '''
    - Plot and save residual of (fit-hist)/hist
    '''
    # Draw residual and get stats of fit
    hres.SetTitle('')
    hres.Draw('e')
    line.SetLineColor(ROOT.kGreen+2)
    ps.c.Update()
    s = hres.GetListOfFunctions().FindObject('stats')
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.75)
    s.SetOptStat(0)
    s.SetOptFit(1111)
    s.Draw()
    # Draw Label
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Draw line at 0
    #l1 = ROOT.TLine(fitlow, 0., fithigh,  0.)
    #l1.SetLineStyle(2)
    #l1.Draw()
    # Save
    ps.save('%s_res'%fitName, log=False,pdf=True)

def plot_fits(fits,int_lumi,titles,extra):
    '''
    - Take in TF1 array of fits; fits = [fit1, fit2, ... , fitN]
      Array of names; names = ['fit1Name','fit2name', ..., 'fitNname']
    - Return a plot of all the fits
    '''
    colors = [ROOT.kOrange+1, ROOT.kBlue, ROOT.kGreen+2, ROOT.kBlack]
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    draw = ''
    for fit,title,color in zip(fits,titles,colors):
        lg.AddEntry(fit,title,'L')
        fit.SetLineColor(color)
        fit.SetTitle('')
        fit.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
        fit.GetXaxis().SetLabelFont(42)
        # assumes original started out with 1 GeV bins, and the xsec is in pb-1. 
        fit.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) 
        fit.Draw(draw)
        draw = 'same'
    lg.Draw()
    ps.save('AllCategoriesFit_%s'%extra,pdf=True,log=True)

def plot_norm_fits(fits,fitlow,fithigh,int_lumi,titles,extra):
    '''
    - Take in TF1 array of fits; fits = [fit1, fit2, ... , fitN]
      Array of names; names = ['fit1Name','fit2name', ..., 'fitNname']
    - Return a plot of all the fits
    '''
    colors = [ROOT.kOrange+1, ROOT.kBlue, ROOT.kGreen+2, ROOT.kBlack]
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    draw = ''
    fitNorms = []
    # Doesn't work unless the list of functions is made before they're drawn?
    for fit in fits:
        # Need to copy the original fit to make everything work
        new = fit.Clone('%s_copy'%fit.GetName())
        fitNorm = make_norm_func(new,fitlow,fithigh)
        fitNorms.append(fitNorm)
    for fit,title,color in zip(fitNorms,titles,colors):
        lg.AddEntry(fit,title,'L')
        fit.SetLineColor(color)
        fit.SetTitle('')
        fit.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
        fit.GetXaxis().SetLabelFont(42)
        fit.GetYaxis().SetTitle('A.U.')
        fit.Draw(draw)
        draw = 'same'
    lg.Draw()
    ps.save('AllCategoriesNormFit_%s'%extra,pdf=True,log=True)

def plot_linear_funcs(uncPlus,uncMinus,uncs,colors,low,fithigh):
    '''
    - Takes in and plots linear functions used for bckg shape robustness test
    '''
    # Linear functions
    lg1 = ROOT.TLegend(0.15,0.65,0.45,0.95)
    lg1.SetFillStyle(0)
    lg1.SetTextFont(42)
    lg1.SetBorderSize(0)
    lg2 = ROOT.TLegend(0.15,0.15,0.45,0.45)
    lg2.SetFillStyle(0)
    lg2.SetTextFont(42)
    lg2.SetBorderSize(0)
    draw = ''
    # Reverse the plus so the legend lines up with the lines nicely
    for plus,uncP,colorP,minus,uncM,colorM in zip(reversed(uncPlus),reversed(uncs),reversed(colors),uncMinus,uncs,colors):
        nameP = '+'+uncP+'%'
        nameM = '-'+uncM+'%'
        lg1.AddEntry(plus,nameP,'L')
        lg2.AddEntry(minus,nameM,'L')
        plus.SetTitle('')
        plus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
        plus.GetXaxis().SetRangeUser(100,fithigh)
        plus.GetYaxis().SetTitle('Scaling')
        plus.SetMinimum(0)
        plus.SetMaximum(2)
        plus.SetLineStyle(1)
        plus.SetLineColor(colorP)
        plus.Draw(draw)
        draw = 'same'
        minus.SetTitle('')
        minus.GetXaxis().SetTitle('m(#mu^{+}#mu^{-}) [GeV]')
        minus.GetXaxis().SetRangeUser(100,fithigh)
        minus.GetYaxis().SetTitle('Scaling')
        minus.SetMinimum(0)
        minus.SetMaximum(2)
        minus.SetLineStyle(2)
        minus.SetLineColor(colorM)
        minus.Draw(draw)
    lg1.Draw()
    lg2.Draw()
    ps.save('linear_tests',log=False,pdf=True)

def plot_syst_func(fit, histMC, uncMinus, uncNames,fitlow, fithigh,cat,title,int_lumi):
    '''
    - Takes in fit, original MC histogram, lists of +/- shape deformation functions
      and names, and some other inputs
    '''
    # Assume uncs is list of 'minus' functions
    for unc,uncName in zip(uncMinus,uncNames):
        plus, minus = make_pm_func(fit,unc,uncName,fitlow,fithigh,cat)
        plot_three(histMC,fit,plus,minus,uncName,cat,title,int_lumi)

#############
# Utilities #
#############

def make_pm_func(fit, unc,uncName,fitlow,fithigh,cat):
    '''
    - Takes in nominal background pdf and uncertainty function 
    - Returns TF1s of plus, minus background pdf shapes
    - Need to test to make naming consistent...
      At the moment the pointer name and the TF1 name need to be
      the same. Ex: func = ROOT.TF1('func',...)
    '''
    # Need to copy the original fit to make everything work
    new = fit.Clone('%s_pm_copy'%fit.GetName())
    plus = ROOT.TF1('plus_%s_%s'%(cat,uncName),'(2-%s)*%s'%(unc.GetName(),new.GetName()),fitlow,fithigh)
    minus = ROOT.TF1('minus_%s_%s'%(cat,uncName),'%s*%s'%(unc.GetName(),new.GetName()),fitlow,fithigh)
    return plus, minus

def make_norm_func(fit,fitlow,fithigh):
    '''
    - Takes in a TF1
    - Returns the same TF1, but normalized to 1
    '''
    NORM = 1/fit.Integral(fitlow,fithigh)
    return ROOT.TF1('%s_Norm'%fit.GetName(),'(%s)*(%s)'%(NORM,fit.GetName()),fitlow,fithigh)

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

def make_k_factor_plot(kFactor,nominal,fitlow,fithigh,int_lumi,cat,title):
    '''
    - Inputs : kFactor function, nominal bckg fit, plot name and settings
    - Outputs : Saved plot of comparison of background shape w and w/o k-factor
    '''
    # Copy nominal fit and set styles
    nomCopy = nominal.Clone('nomCopy_%s'%cat)
    nomCopy.SetTitle('')
    nomCopy.SetLineColor(ROOT.kBlue)
    # New background shape = kfactor function * nominal background shape
    # Set styles
    kFactorBckgShape = ROOT.TF1('kFactorBckgShape','nomCopy_%s*kFactor'%(cat),fitlow,fithigh)
    kFactorBckgShape.SetTitle('')
    kFactorBckgShape.SetLineColor(ROOT.kOrange+1)
    kFactorBckgShape.GetYaxis().SetTitle('Events/%.1f fb^{-1}' % (int_lumi/1000)) 
    kFactorBckgShape.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    # Draw Functions
    kFactorBckgShape.Draw()
    nomCopy.Draw('same')
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(nomCopy,'Nominal','L')
    lg.AddEntry(kFactorBckgShape,'k-factor applied','L')
    lg.Draw()
    # Draw title
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    # Save
    ps.save('%s_kFactor'%cat,log=True)

def plot_k_factor(kFactor,fitlow,fithigh):
    '''
    - Inputs : k-factor function, plot settings
    - Outputs : Saved plot of k-factor function
    '''
    # kfactor function
    # Set styles
    kFactor.SetTitle('')
    kFactor.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    # Draw Function
    kFactor.Draw()
    # Draw Legend
    lg = ROOT.TLegend(0.15,0.15,0.4,0.4)
    lg.SetFillStyle(0)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.AddEntry(kFactor,'k-factor','L')
    lg.Draw()
    # Save
    ps.save('kFactor',log=False,pdf=True)

def plot_ratio(fit1,fit2,fitlow,fithigh,title,cat,extra):
    '''
    - Inputs : 2 fits, title, category
    - Output : plot of fit2/fit1
    '''
    fit1copy = fit1.Clone('%s_copy'%fit1.GetName())
    fit2copy = fit2.Clone('%s_copy'%fit2.GetName())
    ratio = ROOT.TF1('%s_%s_ratio'%(cat,extra),'(%s)/(%s)'%(fit2copy.GetName(),fit1copy.GetName()),fitlow,fithigh)
    #ratio1 = ROOT.TF1('%s_%s_ratio1'%(cat,extra),'%s'%(fit1copy.GetName()),fitlow,fithigh)
    #ratio2 = ROOT.TF1('%s_%s_ratio2'%(cat,extra),'%s'%(fit2copy.GetName()),fitlow,fithigh)
    ratio.SetTitle('')
    ratio.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    ratio.GetYaxis().SetTitle('Ratio of bckg shapes')
    ratio.Draw('l')
    #ratio1.SetTitle('')
    #ratio1.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    #ratio1.GetYaxis().SetTitle('Ratio of bckg shapes')
    #ratio1.SetLineColor(ROOT.kRed)
    #ratio1.Draw('l')
    #ratio2.SetTitle('')
    #ratio2.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) [GeV]')
    #ratio2.GetYaxis().SetTitle('Ratio of bckg shapes')
    #ratio2.SetLineColor(ROOT.kSpring)
    #ratio2.Draw('lsame')
    HeaderLabel = ROOT.TPaveLabel(0.2, 0.92, 0.8, 0.98,title,'NDC')
    HeaderLabel.SetTextFont(42)
    HeaderLabel.SetTextSize(0.8) 
    HeaderLabel.SetBorderSize(0)
    HeaderLabel.SetFillColor(0)
    HeaderLabel.SetFillStyle(0)
    HeaderLabel.Draw()
    ps.save('test_PI_ratio',log=False,pdf=True)
    #ps.save(ratio.GetName(),log=False,pdf=True)

#**************************#
# Main portion of the code #
#**************************#

if __name__=='__main__':
    import sys
    import array as array
    import ROOT
    ROOT.gROOT.SetBatch(True)
    from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
    import argparse

    parser = argparse.ArgumentParser(description='Fit and draw background shape with options for systematic studies')
    parser.add_argument('--tag',dest='tag',help='Tagged directory name for output plots',default='')
    parser.add_argument('--doPI',dest='doPI',action='store_true',help='Make Photon-Induced robustness test plots',default=False)
    parser.add_argument('--doShapeSyst',dest='doShapeSyst',action='store_true',help='Make background shape robustness test plots',default=False)
    parser.add_argument('--doKfactor',dest='doKfactor',action='store_true',help='Factor in k-factor on the background shape',default=False)
    parser.add_argument('--doNonDY',dest='doNonDY',action='store_true',help='Factor in non-DY MC on the background shape',default=True)
    parser.add_argument('--fitData',dest='fitData',action='store_true',help='Fit bckg shape function to data',default=False)
    args = parser.parse_args()
    
    TAG = args.tag
    DOPI = args.doPI
    DOSHAPESYST = args.doShapeSyst
    DOKFACTOR = args.doKfactor
    DONONDY = args.doNonDY
    FITDATA = args.fitData

    set_zp2mu_style()
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ps = plot_saver('plots/'+TAG)
    ROOT.gPad.SetTicks(1)
    ROOT.TH1.AddDirectory(0)

    #hists_dir = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/CMSSW_8_0_3/DataMCSpectraComparison/mc/80X_v1/'
    # Trigger weight applied
    #hists_dir = '/afs/cern.ch/work/f/ferrico/public/Root_ZPrime_DONT_DELETE/DataMCSpectra/mc/'
    # FIXME set hists_dir to where MC files are
    # Trigger weight down applied
    hists_dir = '/afs/cern.ch/work/c/cschnaib/public/mc_trigWeightDown/'
    rebin = 40
    low = 40
    fitlow = 200
    fithigh = 5500
    high = 5500
    int_lumi = 1000.

    fits = []
    if FITDATA: dataFits = []
    titles = []

    if DOSHAPESYST:
        #unc05 = ROOT.TF1('unc05','[0] + [1]*x',low,fithigh)                      
        #unc05.SetParameters(1.0012,-1.71E-5)
        #unc10 = ROOT.TF1('unc10','[0] + [1]*x',low,fithigh)                      
        #unc10.SetParameters(1.0024,-3.41E-5)
        unc15minus = ROOT.TF1('unc15minus','[0] + [1]*x',low,fithigh)
        unc15minus.SetParameters(1.0036,-5.12E-5)
        unc20minus = ROOT.TF1('unc20minus','[0] + [1]*x',low,fithigh)
        unc20minus.SetParameters(1.0048,-6.83E-5)
        unc25minus = ROOT.TF1('unc25minus','[0] + [1]*x',low,fithigh)
        unc25minus.SetParameters(1.0060,-8.53E-5)
        #unc30minus = ROOT.TF1('unc30','[0] + [1]*x',low,fithigh)                      
        #unc30minus.SetParameters(1.0072,-10.24E-5)
        #unc35minus = ROOT.TF1('unc35','[0] + [1]*x',low,fithigh)                      
        #unc35minus.SetParameters(1.0084,-11.94E-5)
        #unc40minus = ROOT.TF1('unc40','[0] + [1]*x',low,fithigh)                      
        #unc40minus.SetParameters(1.0096,-13.65E-5)
        #unc45minus = ROOT.TF1('unc40','[0] + [1]*x',low,fithigh)                      
        #unc45minus.SetParameters(1.0107,-15.36E-5)
        unc50minus = ROOT.TF1('unc50minus','[0] + [1]*x',low,fithigh)
        unc50minus.SetParameters(1.0119,-17.06E-5)
        unc15plus = ROOT.TF1('unc15plus','2-unc15minus',low,fithigh)
        unc20plus = ROOT.TF1('unc20plus','2-unc20minus',low,fithigh)
        unc25plus = ROOT.TF1('unc25plus','2-unc25minus',low,fithigh)
        unc50plus = ROOT.TF1('unc50plus','2-unc50minus',low,fithigh)
        # Lists
        uncPlus =  [ unc15plus, unc20plus, unc25plus, unc50plus]
        uncMinus = [unc15minus,unc20minus,unc25minus,unc50minus]
        uncNames = ['15','20','25','50']
        colors =  [ROOT.kBlack,ROOT.kBlue,ROOT.kGreen+2,ROOT.kOrange+1]
        plot_linear_funcs(uncPlus,uncMinus,uncNames,colors,low,fithigh)

    category = {'BB':('Our2012MuonsPlusMuonsMinusBarrelHistos','DimuonMassVertexConstrained_bb','Barrel - Barrel'),
                'BEp':('Our2012MuonsPlusMuonsMinusNegativeHistos','DimuonMassVertexConstrained_ne','Barrel - Positive Endcap'),
                'BEn':('Our2012MuonsPlusMuonsMinusPositiveHistos','DimuonMassVertexConstrained_pe','Barrel - Negative Endcap'),
                'All':('Our2012MuonsPlusMuonsMinusHistos','DimuonMassVertexConstrained','All Categories')
               }
    toDo = ['BB','BEp','BEn','All']

    for cat in toDo:
        rootDir,rootHist,title = category[cat]
        print rootDir, rootHist

        # Make nominal MC histogram
        histMC = make_mc_hist(DONONDY,False,int_lumi,low,high,rebin,hists_dir,rootDir,rootHist,title,cat)
        # Make nominal MC fit and saves plot
        nominal = fit_bckg(histMC,cat+'_Fit',low,high,fitlow,fithigh,rebin,title)
        fits.append(nominal)
        titles.append(title)

        # Include Photon-Induced background
        if DOPI:
            # Make MC histogram with PI background
            histMCPI = make_mc_hist(DONONDY,DOPI,int_lumi,low,high,rebin,hists_dir,rootDir,rootHist,title,cat)
            nominalPI = fit_bckg_PI(histMCPI,cat+'_PI_Fit',low,high,fitlow,fithigh,rebin,title)
            # Plot Bckg and Bckg+PI fit functions
            plot_two(nominal,'Nominal',histMC,nominalPI,'Nominal+PI',histMCPI,'nominal_PI','Fit',title,cat,int_lumi)
            # Plot Ratio Bckg+PI / Bckg
            plot_ratio(nominal,nominalPI,fitlow,fithigh,title,cat,'PI')

        # Plot +/nominal/- background shapes with linear deformation functions
        if DOSHAPESYST: plot_syst_func(nominal,histMC,uncMinus,uncNames,fitlow,fithigh,cat,title,int_lumi)
        
        # Factor in k-factor parametrization from Sam
        if DOKFACTOR: 
            kFactor = ROOT.TF1('kFactor','TMath::Min(1.01696 - 7.73522E-5*x + 6.69239E-9*x*x,1.0)',low,high)
            make_k_factor_plot(kFactor,nominal,fitlow,fithigh,int_lumi,cat,title)

        if FITDATA: 
            data_dir = '/afs/cern.ch/work/c/cschnaib/public/'
            dataRebin = 50
            dataFitHigh = 1500
            dataHigh = 2000
            dataHist = make_data_hist(data_dir,rootDir,rootHist,dataRebin,low,dataHigh,title,cat,13000)#,varBins=True)
            dataFit = fit_data(dataHist,'%s_dataFit'%(cat),low,dataHigh,fitlow,dataFitHigh,dataRebin,title)
            dataFits.append(dataFit)

        print '\n*****************\n'
    print '\n'

    # Plot All categories on top of each other
    plot_fits(fits,int_lumi,titles,'mc')
    plot_norm_fits(fits,fitlow,fithigh,int_lumi,titles,'mc')
    if FITDATA:
        plot_fits(dataFits,13000,titles,'data')
        plot_norm_fits(dataFits,fitlow,1500,13000,titles,'data')
    # Plot k-factor function
    if DOKFACTOR: plot_k_factor(kFactor,fitlow,fithigh)
    print 'done'
    


