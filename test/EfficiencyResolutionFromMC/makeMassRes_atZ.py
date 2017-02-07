#!/usr/bin/python

# import ROOT in batch mode
import sys,os
import argparse
import math

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ptbins = [32, 52, 72, 100, 152, 200, 300, 452, 800]

ROOT.gROOT.LoadMacro("cruijff.C+")

def getBinRange(histo,mlow,mhigh):
    ymin =  0
    ymax = -1
    nbins = histo.GetNbinsY()
    for bin in range(nbins):
        if mlow==histo.GetYaxis().GetBinLowEdge(bin): 
            ymin = bin
        if mhigh==histo.GetYaxis().GetBinLowEdge(bin):
            ymax = bin
    return ymin,ymax

def loadHistos(inputdata,inputMC,region,rebin=1):
    _fileDATA = ROOT.TFile(inputdata)
    _fileMC   = ROOT.TFile(inputMC)
    
    hdata = ROOT.TH2F()
    hmc   = ROOT.TH2F()
    
    hdata.SetDirectory(0)
    hmc  .SetDirectory(0)
    ROOT.TH1.AddDirectory(ROOT.kFALSE)
    
    reg = ""
    if   ("BB" in region or "barrel" in region):  reg = "_BB"
    elif ("BE" in region or "endcap" in region):  reg = "_BE"

    hdata = _fileDATA.Get("Resolutiontunepnew/DileptonMass_2d_vsPt%s" %(reg)).Clone()
    hmc   = _fileMC  .Get("Resolutiontunepnew/DileptonMass_2d_vsPt%s" %(reg)).Clone()

    _fileDATA.Close()
    _fileMC  .Close()

    data = [ROOT.TH1D() for x in range(len(ptbins)-1)]
    mc   = [ROOT.TH1D() for x in range(len(ptbins)-1)]
    ptda = [0 for x in range(len(ptbins)-1)]
    ptmc = [0 for x in range(len(ptbins)-1)]
    for h in data:
        h.SetDirectory(0)
        ROOT.TH1.AddDirectory(ROOT.kFALSE)
    for h in mc:
        h.SetDirectory(0)
        ROOT.TH1.AddDirectory(ROOT.kFALSE)
    
    for i,h in enumerate(data):
        ymin,ymax=getBinRange(hdata,ptbins[i],ptbins[i+1])
        hdata.GetYaxis().SetRange(ymin,ymax)
        hmc.GetYaxis().SetRange(ymin,ymax)
        ptda[i] = hdata.GetMean(2)
        ptmc[i] = hmc.GetMean(2)
        hdata.GetYaxis().SetRange()
        hmc  .GetYaxis().SetRange()

        data[i] = hdata.ProjectionX("datapy%s%s" %(ptbins[i],region),ymin,ymax)
        mc  [i] = hmc  .ProjectionX("mcpy%s%s" %(ptbins[i],region)  ,ymin,ymax)
        
        data[i].Rebin(2)
        mc  [i].Rebin(2)
        
        if (data[i].Integral() < 1000): 
            data[i].Rebin(2)
            mc[i].Rebin(2)
        
        if (mc[i].Integral() < 200):
            mc[i].Rebin(2)

    return data,mc,ptda,ptmc

def doFit(hist,output,rap="BB",flavour="DATA"):
    c1 = ROOT.TCanvas("c1","c1",700,700)
    c1.cd()

    sig    = []
    sige   = []
    mean   = []
    meane  = []
    
    DOCRYSTALBALL = False
    DOCRUIJFF = True
    
    for i,h in enumerate(hist):
        print "+++++++++++++++++++++++++++++++++++++++++"
        print "Fitting histogram for %d < pt_{l} <%d" %(ptbins[i],ptbins[i+1])
        print "+++++++++++++++++++++++++++++++++++++++++\n"

        fit_min = h.GetMean() - h.GetRMS() 
        fit_max = h.GetMean() + h.GetRMS()

        # fit with a gaussian 
        gaus = ROOT.TF1("gaus","gaus",fit_min,fit_max)
        gaus.SetLineColor(ROOT.kBlue)
        gaus.SetParameters(0,h.GetMean(),h.GetRMS())
        h.Fit("gaus","M0R+")

        if DOCRYSTALBALL:
            crystal = ROOT.TF1("crystal","crystalball",fit_min,fit_max)
            crystal.SetLineColor(ROOT.kRed)
        
            tmp_mean = gaus.GetParameter(1)
            tmp_sig  = gaus.GetParameter(2)
            min_mean = tmp_mean+tmp_mean
            max_mean = tmp_mean-tmp_mean
            crystal.SetParameters(gaus.GetParameter(0), tmp_mean, tmp_sig, 1.4, 2)
            crystal.SetParLimits(0, 0, 2*gaus.GetParameter(0)) # //const
            crystal.SetParLimits(1, min_mean, max_mean)
            crystal.SetParLimits(2, 0, 2.5*h.GetRMS())# //sigma
            crystal.SetParLimits(3, 0.5, 2.)# //alpha
            crystal.SetParLimits(4, 0., 3.)# //alpha            
            h.Fit("crystal","M0R+")
        
            # crystal
            mean  .append(crystal.GetParameter(1))
            meane .append(crystal.GetParError(1))
            sig   .append(crystal.GetParameter(2))
            sige  .append(crystal.GetParError(2))
        
        elif DOCRUIJFF:
            fit_min = 60.
            fit_max = 120. 
            if ptbins[i] < 100: fit_max = 105.

            cruijff = ROOT.TF1("cruijff",ROOT.cruijff,fit_min,fit_max,5)
            cruijff.SetParNames("Constant","Mean","Sigma","AlfaL","AlfaR")
            cruijff.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1),gaus.GetParameter(2),1.4,1.3)
            cruijff.SetLineColor(ROOT.kRed)
            h.Fit("cruijff","M0R+")
            
            mean  .append(cruijff.GetParameter(1))
            meane .append(cruijff.GetParError(1))
            sig   .append(cruijff.GetParameter(2))
            sige  .append(cruijff.GetParError(2))
            
        else:
            mean  .append(gaus.GetParameter(1))
            meane .append(gaus.GetParError(1))
            sig   .append(gaus.GetParameter(2))
            sige  .append(gaus.GetParError(2))
  
#        alpha .append(crystal.GetParameter(3))

#        alphae.append(crystal.GetParError(3))
#        n     .append(crystal.GetParameter(4))  #exponential parameter 
#        ne    .append(crystal.GetParError(4))   #exponential parameter

        h.SetTitle("Mass resolution for %d < p_{T} <%d" %(ptbins[i],ptbins[i+1]))
        h.GetXaxis().SetTitle("m_{ll} [GeV]")
        h.SetLineColor(ROOT.kBlack)
        h.GetXaxis().SetRangeUser(fit_min,fit_max)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.7)
                
        h.Draw("E")
        if DOCRYSTALBALL:
            crystal.Draw("SAME")
        elif DOCRUIJFF:
            cruijff.Draw("SAME")
        else:
            gaus.Draw("SAME")
            
        saveas = "/MassRes_%s_Pt%d_%d_%s" %(flavour,ptbins[i],ptbins[i+1],rap)
        c1.SaveAs(output+saveas+".root")
        c1.SaveAs(output+saveas+".C")
        c1.SaveAs(output+saveas+".png")
        c1.SaveAs(output+saveas+".pdf")

    print "DONE Fitting..."
    return mean,meane,sig,sige

    
def drawMassRes(data,mc,output,rapidity,ptda,ptmc):
    
    pt_e = [0 for x in range(len(ptbins)-1)]
    pt_x = [0 for x in range(len(ptbins)-1)]
    for i,pt in enumerate(pt_x):
        pt_x[i] = ptbins[i]+(ptbins[i+1]-ptbins[i])/2.
        pt_e[i] = (ptbins[i+1]-ptbins[i])/2.

    (da_mean,da_meane,da_sig,da_sige) = doFit(data,output,rapidity,"DATA")
    (mc_mean,mc_meane,mc_sig,mc_sige) = doFit(mc  ,output,rapidity,"MC")
    
    c2 = ROOT.TCanvas("c2","c2",700,700)
    c2.cd()

    # Upper plot will be in pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetGrid()        # Vertical grid
    pad1.SetBottomMargin(0)
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()               # pad1 becomes the current pad
    pad1.SetTicks()
    #    fun = ROOT.TF1("fun","pol2")
    #    fun.SetParameters(0.,0.,0.)
  
    res_data  = ROOT.TGraphAsymmErrors(len(pt_x))
    res_data.SetName("res_data")
    res_mc    = ROOT.TGraphAsymmErrors(len(pt_x))
    res_mc  .SetName("res_mc")
    ratio     = ROOT.TGraphErrors(len(pt_x))
    ratio   .SetName("ratio")
    for i,pt in enumerate(pt_x):        
        res_data.SetPoint(i,ptda[i],da_sig[i])
        res_data.SetPointError(i,ptda[i]-ptbins[i],ptbins[i+1]-ptda[i],da_sige[i],da_sige[i])
        res_mc  .SetPoint(i,ptmc[i],mc_sig[i])
        res_mc  .SetPointError(i,ptmc[i]-ptbins[i],ptbins[i+1]-ptmc[i],mc_sige[i],mc_sige[i])
        ratio   .SetPoint(i,pt,da_sig[i]/mc_sig[i])
        ratio   .SetPointError(i,pt_e[i],(da_sig[i]/mc_sig[i])*math.sqrt((da_sige[i]/da_sig[i])**2+(mc_sige[i]/mc_sig[i])**2))
        
    res_data.SetMarkerStyle(22)
    res_data.SetMarkerColor(ROOT.kBlack)
    res_data.SetLineColor(ROOT.kBlack)
    res_data.SetFillColor(0)
    res_data.SetTitle("Dimuon mass resolution vs pT")
    res_data.GetYaxis().SetTitle("#sigma (Z peak)")
    res_data.GetYaxis().SetTitleOffset(1.2)
#    res_data.GetXaxis().SetTitle("p_{T} (#mu^{#pm}) [GeV]")
    res_data.GetYaxis().SetRangeUser(0.1,8)
    res_data.GetXaxis().SetRangeUser(30,800)
#    fun.SetParameters(0.,0.,0.)
#    res_crystal.Fit(fun,"M+")
#    res_crystal.GetFunction("fun").SetLineColor(ROOT.kRed+2)
    res_data.Draw("AP E0")
    
    

    res_mc.SetMarkerStyle(22)
    res_mc.SetMarkerColor(ROOT.kRed)
    res_mc.SetLineColor(ROOT.kRed)
    res_mc.SetFillColor(0)
    res_mc.SetTitle("Dimuon mass resolution vs pT")
    res_mc.GetYaxis().SetTitle("Resolution (Z peak)")
    res_mc.GetYaxis().SetTitleOffset(1.5)
#    res_mc.GetXaxis().SetTitle("p_{T} (#mu^{#pm}) [GeV]")
#    res_mc.GetYaxis().SetRangeUser(0,.1)
#    res_mc.GetXaxis().SetRangeUser(0,600)
#    fun.SetParameters(0.,0.,0.)
#    res_crystal.Fit(fun,"M+")
#    res_crystal.GetFunction("fun").SetLineColor(ROOT.kRed+2)
    res_mc.Draw("P E0 SAME")

    leg = ROOT.TLegend(.35,.7,.50,.80,"","brNDC")
    leg.AddEntry(res_data,"DATA")
    leg.AddEntry(res_mc,"Simulation")
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.SetTextSize(.02)
    leg.Draw("SAME")

    c2.cd()          # Go back to the main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2",0, 0.1, 1, 0.30)    
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.25)
    pad2.SetGrid()
    pad2.Draw()
    pad2.cd()
    pad2.SetTicks()

    ratio.SetMarkerColor(ROOT.kBlue-4)
    ratio.SetFillColor(ROOT.kBlue-4 )
    ratio.SetTitle("")
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetXaxis().SetTitle("p_{T} (#mu^{#pm}) [GeV]")
    ratio.GetYaxis().SetRangeUser(0.5,1.5)
    ratio.GetXaxis().SetRangeUser(30,800)
    ratio.GetYaxis().SetTitleOffset(0.40)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetLabelSize(0.14)    
    ratio.GetYaxis().SetNdivisions(506)    
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetLabelSize(0.14)
            
    ratio.Draw("A P E2")
    pad2.Update()
    line = ROOT.TLine(30,1,800,1)

    line.SetLineColor(ROOT.kBlue+1)
    line.SetLineWidth(2)
    line.Draw()

    saveas = "/MassResolutionVsPt_%s" %(rapidity)
    c2.SaveAs(output+saveas+".png")
    c2.SaveAs(output+saveas+".pdf")
    c2.SaveAs(output+saveas+".root")
    c2.SaveAs(output+saveas+".C")
    
    # PRINT FIT RESULTS!!!
#    ndf = []
#    chi = []
#    for h in hist:
#        chi.append(h.GetFunction("crystal_%s"%nrms).GetChisquare())
#        ndf.append(h.GetFunction("crystal_%s"%nrms).GetNDF())

    print "|--------------------------------------------------------------------------|"
    print "|                   MASS RESOLUTION PARAMETRIZATION                        |"
    print "|--------------------------------------------------------------------------|"
    print "|      pT   %s       |      Mean:  Data / MC      |    Sigma: Data / MC     | " %(rapidity)
    print "|--------------------------------------------------------------------------|"
    
    for i in range(0,len(ptbins)-1):
        if ptbins[i]>200 and rapidity=="E": break
        print "| %4d < p_T < %4d |  %5.3f / %5.3f ( %5.3f ) | %5.3f / %5.3f ( %5.3f ) |" %(pt_x[i]-pt_e[i], pt_x[i]+pt_e[i],
	                                                                                   da_mean[i], mc_mean[i], da_mean[i]/mc_mean[i],
                                                                                           da_sig[i], mc_sig[i], da_sig[i]/mc_sig[i])
    print "|--------------------------------------------------------------------------|"

    
def makeMassRes(inputDATA,inputMC,output):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetStatX(.9)
    ROOT.gStyle.SetStatY(.9)
    
    (data_B,mc_B,ptda,ptmc) = loadHistos(inputDATA,inputMC,"BB")
 #   (data_E,mc_E) = loadHistos(inputDATA,inputMC,"BE")

    drawMassRes(data_B,mc_B,output,"BB",ptda,ptmc)
#    drawMassRes(data_E,mc_E,output,"BE")
    
#    
#    resBB.SetMarkerStyle(22)
#    resBB.SetMarkerColor(ROOT.kRed)
#    resBB.SetLineColor(ROOT.kRed)
#    resBB.SetFillColor(0)
#    resBB.SetTitle("Dimuon mass resolution vs pT")
#    resBB.GetYaxis().SetTitle("Dimuon Mass Resolution")
#    resBB.GetYaxis().SetTitleOffset(1.5)
#    resBB.GetXaxis().SetTitle("p_T (#mu^{#pm}) [GeV]")
#    resBB.GetYaxis().SetRangeUser(0,.2)
#    resBB.GetXaxis().SetRangeUser(0,5000)
#    resBB.GetFunction("fun").SetLineColor(ROOT.kRed+1)
#    resBB.Draw("AP E0")
#    
#    resBE.SetMarkerStyle(22)
#    resBE.SetMarkerColor(ROOT.kBlue+1)
#    resBE.SetLineColor(ROOT.kBlue+1)
#    resBE.SetFillColor(0)
#    resBE.SetTitle("Dimuon mass resolution vs mass")
#    resBE.GetYaxis().SetTitle("Dimuon Mass Resolution")
#    resBE.GetYaxis().SetTitleOffset(1.5)
#    resBE.GetXaxis().SetTitle("p_T (#mu^{#pm}) [GeV]")
#    resBE.GetYaxis().SetRangeUser(0,.2)
#    resBE.GetXaxis().SetRangeUser(0,5000)
#    resBE.GetFunction("fun").SetLineColor(ROOT.kBlue+2)
#    resBE.Draw("PE0 SAME")
#        
#    leg = ROOT.TLegend(.35,.7,.50,.80,"","brNDC")
#    leg.AddEntry(resBB,"BB")
#    leg.AddEntry(resBE,"BE")
#    leg.SetTextFont(42)
#    leg.SetBorderSize(0)
#    leg.SetTextSize(.02)
#    leg.Draw("SAME")
#    
#    res.SetGrid()
#    saveas = "/MassResolutionVsPt_2CAT"
#    res.SaveAs(output+saveas+".png")
#    res.SaveAs(output+saveas+".pdf")
#    res.SaveAs(output+saveas+".root")
#    res.SaveAs(output+saveas+".C")
    
         
#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="makeMassRes.py [options]",description="Compute mass resolution",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--iDATA", dest="inputDATA",default="", help='Input filename')
    parser.add_argument("--iMC", dest="inputMC",default="", help='Input filename')
    parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name to store results')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output);
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+args.output)
    
    inputDATA = args.inputDATA
    inputMC   = args.inputMC
    output=args.output
    
    print "Running on: %s %s with 2 categories" %(inputDATA,inputMC)
    print "Saving result in: %s" %(output)

    makeMassRes(inputDATA,inputMC,output)
    print "DONE"
