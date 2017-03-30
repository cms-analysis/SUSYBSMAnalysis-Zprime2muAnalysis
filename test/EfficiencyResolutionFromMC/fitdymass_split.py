#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.TH1.AddDirectory(0)


# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(111)

# variable = 'DimuonMassVertexConstrained'
variable = 'DimuonMassVertexConstrained_bb'
# variable = 'DimuonMassVertexConstrained_be'

low = fitlow = 150
high = fithigh = 5000
high_for_data = 2500

# ps = plot_saver('plots/fitdymass/split'+ variable)
ps = plot_saver('plots/'+ variable)


int_lumi = 36295.39
rebin = 40
use_non_dy = True

# Masses = [50, 120, 200, 400, 800, 1400, 2300, 3500, 4500]


masses  = ['dy50to120', 'dy120to200', 'dy200to400', 'dy400to800', 'dy800to1400', 'dy1400to2300', 'dy2300to3500', 'dy3500to4500', 'dy4500to6000']
nevents = [2977600, 100000, 100000,  98400, 100000, 100000, 100000, 100000, 100000]
#sigmas  = [  1915.,  12.2,  1.53, 0.0462, 0.00586, 0.00194, 1.70e-4, 2.21e-5] # in pb, PYTHIA*1.3
sigmas  = [  1975, 19.32, 2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7]
weights = [int_lumi / nev * sig for nev,sig in zip(nevents, sigmas)]
#weights = [x/weights[-1] for x in weights]

hists = []
hists_dir = '../DataMCSpectraComparison/mc/'
for m,w in zip(masses, weights):
    fn = 'ana_datamc_%s.root' % m #if m != 20 else 'ana_datamc_zmumu.root'
    fn = hists_dir + fn
    f = ROOT.TFile(fn)
    d = f.Our2016MuonsPlusMuonsMinusHistos
    h = d.Get(variable).Clone('%s' % m)
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(low, high)
    print m,w
    h.Scale(w)
    h.Draw()
    ps.save('rawmass%s' % m)
    hists.append(h)
 
if use_non_dy:
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
    non_dy_samples = [WWinclusive, WW200to600, WW600to1200, WW1200to2500, WW2500, 
    						WZ_skim,
							ZZ_skim,
    						WZ_ext, 
    						ZZ_ext_skim,
    						Wantitop, tW, 
    						Wjets, 
    						ttbar_lep50to500, 
							ttbar_lep_500to800, 
							ttbar_lep_800to1200, 
							ttbar_lep_1200to1800, 
							ttbar_lep1800toInf,
							#     						qcd50to80, 
    						qcd80to120, qcd120to170, 
    						qcd170to300, 
    						qcd300to470, qcd470to600, qcd600to800, qcd800to1000, qcd1000to1400, 
    						qcd1400to1800, qcd1800to2400, qcd2400to3200, qcd3200, 
    						dyInclusive50
    						]
    for sample in non_dy_samples:
        fn = 'ana_datamc_%s.root' % sample.name
        fn = hists_dir + fn
        f = ROOT.TFile(fn)
        d = f.Our2016MuonsPlusMuonsMinusHistos

        w = sample.partial_weight * int_lumi
        h = d.Get(variable).Clone('%s' % sample.name)
        h.Rebin(rebin)
        h.GetXaxis().SetRangeUser(low, high)
        h.Scale(w)
        h.Draw()
        ps.save('rawmass_%s' % sample.name)
        print sample.name, w
        hists.append(h)
  
  
htot = hists[0].Clone()

for j in xrange(1, len(hists)):
    htot.Add(hists[j])

htot.SetTitle('')
htot.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
htot.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.

htot.Draw()
    
def fit_it(lo, hi):
    
    print " ----------------------------------------------------------------------------------====================== INIZIO "

    htot.Draw()
    ps.c.Update()

#     fcn = ROOT.TF1("fcn", "(x<=500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
#     						(x>500)*(exp([4] + [5]*x + [6]*x*x)*x**([7]))", lo, hi)
#     fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "kH")
#     fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -3.5)
    fcn = ROOT.TF1("fcn", "(x <= 500)*(exp([0] + [1]*x + [2]*x*x)*x**([3])) + \
    					   (x> 500)*(exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)*x**([8]))", lo, hi)
    fcn.SetParNames("aL", "bL", "cL", "kL", "aH", "bH", "cH", "dH", "kH")
    fcn.SetParameters(24, -2E-3, -1E-7, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)



#     fcn = ROOT.TF1("fcn", "(x<=500)*(exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])) + \
#     						(x>500)*(exp([5] + [6]*x + [7]*x*x + [8]*x*x*x)*x**([9]))", lo, hi)
#     fcn.SetParNames("aL", "bL", "cL", "dL", "kL", "aH", "bH", "cH", "dH", "kH")
#     fcn.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5, 24, -2E-3, -1E-7, -1E-10, -3.5)

#     fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", lo, hi)
#     fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", 500, hi)
#     fcn = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", lo, 500)
#     fcn.SetParNames("a", "b", "c", "d", "k")
#     fcn.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5)
    fcn.SetLineColor(ROOT.kBlue)
    htot.Fit(fcn, 'SREMV')    
    s = htot.GetListOfFunctions().FindObject("stats")
    s.SetName("Const")
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.69)
    s.SetY2NDC(0.94)
    s.SetOptStat(10)
    s.SetOptFit(11111)
    s.SetTextColor(ROOT.kBlue)


#     fcn_2 = ROOT.TF1("fcn", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)*x**([4])", 500, hi)
#     fcn_2.SetParNames("a", "b", "c", "d", "k")
#     fcn_2.SetParameters(24, -2E-3, -1E-7, -1E-10, -3.5)
#     fcn_2.SetLineColor(ROOT.kRed)
#     htot.Fit(fcn_2, 'SREMV')
#     htot.Fit(fcn, 'SREMV')  
#     s_2 = htot.GetListOfFunctions().FindObject("stats")
#     s_2.SetName("Const")
#     s_2.SetX1NDC(0.73)
#     s_2.SetY1NDC(0.44)
#     s_2.SetY2NDC(0.69)
#     s_2.SetOptStat(10)
#     s_2.SetOptFit(11111)
#     s_2.SetTextColor(ROOT.kRed)


    ps.c.Update()   

#     s.Draw("e")
#     s_2.Draw("e")

    ps.save('mass%i_%i' % (lo, hi), pdf=True, pdf_log=True)
    
    htot.GetXaxis().SetRangeUser(0,500)
    ps.save('mass%i_%i_ZOOM_500' % (lo, hi), pdf=True, pdf_log=True)
    
    htot.GetXaxis().SetRangeUser(0,1000)
    ps.save('mass%i_%i_ZOOM_1000' % (lo, hi), pdf=True, pdf_log=True)
        
    ps.c.Update()
	

    xax = htot.GetXaxis()
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-hist)/hist', htot.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(htot.GetNbinsX()+1))
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, hres.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= lo and xhi <= hi:#and xhi < 500:
            res = fcn.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
            if htot.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/htot.GetBinContent(i))
                hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))
#         if xlo >= lo and xhi <= hi and xhi > 500:
#             res = fcn_2.Integral(xlo, xhi)/(xhi-xlo) - htot.GetBinContent(i)
#             if htot.GetBinContent(i) > 0:
#                 hres.SetBinContent(i, res/htot.GetBinContent(i))
#                 hres.SetBinError(i, htot.GetBinError(i)/htot.GetBinContent(i))
     
    hres.SetMinimum(-0.25)
    hres.SetMaximum(0.25)
    hres.GetXaxis().SetRangeUser(lo, hi)
    hres.Draw('e')
    l1 = ROOT.TLine(lo, 0., hi,  0.)
    l1.Draw()
    
    t = ROOT.TPaveLabel(0.15, 0.825, 0.45, 0.925, "K factor: new", 'brNDC')
    t.SetTextFont(42)
    t.SetTextSize(0.5)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetFillStyle(0)
#     t.Draw()
    
    ps.save('res_%i_%i' % (lo, hi), pdf=True, log=False)
    
    ps.c.Update()
       
       
    f_data = ROOT.TFile('/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/data/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root')
    c_data = f_data.Our2016MuonsPlusMuonsMinusHistos
    data = c_data.Get(variable)
    data.Rebin(rebin)
    data.SetStats(0)
    data.GetXaxis().SetTitle('reconstructed m(#mu^{+}#mu^{-}) (GeV)')
    data.GetYaxis().SetTitle('Events/%i GeV/%.1f fb^{-1}' % (rebin, int_lumi/1000)) # assumes original started out with 1 GeV bins, and the xsec is in pb-1.
    data.Draw()
    htot.SetLineColor(ROOT.kRed)
    htot.GetXaxis().SetRangeUser(lo, hi)
    htot.Draw("same")
    htot.Fit(fcn, 'SREMV')
    htot.GetXaxis().SetRangeUser(lo, high_for_data)
    data.GetXaxis().SetRangeUser(lo, high_for_data)
#     htot.SetMinimum(10e-5)
#     htot.SetMaximum(10e5)
    data.SetMinimum(10e-6)
    data.SetMaximum(10e5)
    ps.save('data_%i_%i' % (lo, hi), pdf_log=True)
    
    data.SetMinimum(10)
    data.GetXaxis().SetRangeUser(0,500)
    ps.save('data%i_%i_ZOOM_500' % (lo, hi), pdf_log=True)

    data.SetMinimum(1)    
    data.GetXaxis().SetRangeUser(0,1000)
    ps.save('data%i_%i_ZOOM_1000' % (lo, hi), pdf_log=True)
    
    
    
    data.GetXaxis().SetRangeUser(lo, high_for_data)
    hres = ROOT.TH1F('hres_%i_%i' % (lo, hi), ';m(#mu^{+}#mu^{-}) (GeV);(fit-data)/data', data.GetNbinsX(), xax.GetBinLowEdge(1), xax.GetBinLowEdge(data.GetNbinsX()+1))
    xax = data.GetXaxis()
    for h in [hres]:
#        h.GetYaxis().SetLabelSize(0.02)
        h.SetMarkerStyle(2)
        h.SetStats(0)
    for i in xrange(1, data.GetNbinsX()+1):
        xlo = xax.GetBinLowEdge(i)
        xhi = xax.GetBinLowEdge(i+1)
        if xlo >= lo and xhi <= hi:
#             print data.GetBinContent(i)
            res = fcn.Integral(xlo, xhi)/(xhi-xlo) - data.GetBinContent(i)
#             print res, data.GetBinContent(i), xlo, xhi
            if data.GetBinContent(i) > 0:
                hres.SetBinContent(i, res/data.GetBinContent(i))
                hres.SetBinError(i, data.GetBinError(i)/data.GetBinContent(i))
            if data.GetBinContent(i) == 0:
                hres.SetBinContent(i, 1)
                hres.SetBinError(i, 0)


    hres.SetMinimum(-1)
    hres.SetMaximum( 1)
    hres.GetXaxis().SetRangeUser(lo, high_for_data)
    hres.Draw('e')
    l1 = ROOT.TLine(lo, 0., high_for_data,  0.)
    l1.Draw()

#     t = ROOT.TPaveLabel(0.15, 0.825, 0.45, 0.925, "K factor: old", 'brNDC')
#     t.SetTextFont(42)
#     t.SetTextSize(0.5)
#     t.SetBorderSize(0)
#     t.SetFillColor(0)
#     t.SetFillStyle(0)
#     t.Draw()
    
    ps.save('res_%i_%i_data' % (lo, hi), pdf=True, log=False)
 
    ps.c.Update() 
 
 
 
    
#    print fcn.GetProb(), fcn.GetParameter(2), fcn.GetParError(2)
#     c = r.GetCovarianceMatrix()
#     r.Print("V")
#     print " ----------------------------------------------------------------------------------====================== FINE "
#     print c
	
	

#l = range(200, 2200, 200)
# l = [60, 120, 140, 160, 200]
l = [150]
for lo in l:
    print " ----------------------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>> %d " %lo
    fit_it(lo, high)
    
print " ------------------------------------------------------------------------- Passo al Draw "
print variable, high
#fit_it(400,2000)

# Take different Drell-Yan mass spectra and plot them overlayed,
# then calculate and plot their ratios.
#fsets = ["", "_c1", "_c2"]
# draw_overlay(fsets)
