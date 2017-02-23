#!/usr/bin/python

# import ROOT in batch mode
import sys,os
import argparse
import math
from setTDRStyle import setTDRStyle

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv
    
mrange = [120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500]
#mrange = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500]

ROOT.gROOT.LoadMacro("cruijff.C+")

FITMIN = -1.
FITMAX =  1. 
def getBinRange(histo,mlow,mhigh):
    xmin =  0
    xmax = -1
    nbins = histo.GetNbinsX()
    for bin in range(nbins):
        if mlow==histo.GetXaxis().GetBinLowEdge(bin): 
            xmin = bin
        if mhigh==histo.GetXaxis().GetBinLowEdge(bin):
            xmax = bin
    return xmin,xmax
    
def loadHistos(inputfile,region,rebin):
    _file = ROOT.TFile(inputfile)

    res = ROOT.TH2F()
    res.SetDirectory(0)
    ROOT.TH1.AddDirectory(ROOT.kFALSE)    
    histoname = "Resolutiontunepnew"
    if ("BB" in region):
        res = _file.Get("%s/DileptonMassResVMass_2d_BB" %(histoname)).Clone()
    elif ("BE" in region):
        res = _file.Get("%s/DileptonMassResVMass_2d_BE" %(histoname)).Clone()

    histos = [ROOT.TH1D() for x in range(len(mrange)-1)]
    for h in histos:
        h.SetDirectory(0)
        ROOT.TH1.AddDirectory(ROOT.kFALSE)    
        
    c1 = ROOT.TCanvas("c1","c1",700,700)
    c1.cd()

    for i,h in enumerate(histos): 
        xmin,xmax = getBinRange(res,mrange[i],mrange[i+1])
        histos[i] = res.ProjectionY("res%s%s" %(mrange[i],region), xmin, xmax)
        histos[i].Rebin(5)
        
#        if (histos[i].Integral() < 5000): 
#            histos[i].Rebin(2)
        

        
    return histos

def doFitGeneric(hist,output,rap="BB",fit="cruijff",syst=False):
    c1 = ROOT.TCanvas("c1","c1",700,700)
    c1.cd()

    pars = []
    errs = []
    chi2 = []
    for i,h in enumerate(hist):
        if ("cruijff" in fit):
            fit_min = h.GetMean()-2.0*h.GetRMS()
            fit_max = h.GetMean()+1.7*h.GetRMS()
        elif ("crystal" in fit):
            fit_min = h.GetMean() - 2.3*h.GetRMS() 
            fit_max = h.GetMean() + 1.0*h.GetRMS()
        else: 
            fit_min = FITMIN
            fit_max = FITMAX

        print "+++++++++++++++++++++++++++++++++++++++++"
        print "Fitting histogram for %d < m_{ll} <%d, with Range=[%3.2f, %3.2f]" %(mrange[i],mrange[i+1],fit_min,fit_max)
        print "+++++++++++++++++++++++++++++++++++++++++\n"

 
        # fit with a gaussian to use parameters of the fit for the CB...
        gaus = ROOT.TF1("gaus","gaus",fit_min,fit_max)
        gaus.SetParameters(0,h.GetMean(),h.GetRMS())
        h.Fit("gaus","M0R+")
        
        funct = ROOT.TF1()
        if "cruijff" in fit: 
            print ">>>>>> Using Cruijff >>>>>>>>"
            funct = ROOT.TF1(fit,ROOT.cruijff,fit_min,fit_max,5)
            funct.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1), gaus.GetParameter(2), 0., 0.) #15, 0.001)             
            funct.SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR")        

        elif "crystal" in fit: 
            print ">>>>>>>>  Using CRYSTAL BALL >>>>>>>>"
            funct = ROOT.TF1(fit,"crystalball",fit_min,fit_max)
            funct.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1), gaus.GetParameter(2), 1.4, 1.5)
#            funct.SetParLimits(1, gaus.GetParameter(1)*0.5, gaus.GetParameter(1)*1.5)
            funct.SetParLimits(2, 0., 2.5*h.GetRMS())
            funct.SetParLimits(3, 0.5, 2.)
            funct.SetParLimits(4, 0., 3.)
            
        funct.SetLineColor(ROOT.kBlue)
        funct.SetLineWidth(2)
        h.Fit(fit,"M0R+")

        if syst: 
            if ("cruijff" in fit): 
                print ">>>>>>>>  Using CRYSTAL BALL for systematics >>>>>>>>"
                systfunc = ROOT.TF1("systfunc","crystalball",h.GetMean() - 2.5*h.GetRMS(),h.GetMean() + 1.2*h.GetRMS())
                systfunc.SetParameters(funct.GetParameter(0), funct.GetParameter(1), funct.GetParameter(2), 1.4, 2.)
            elif ("crystal" in fit): 
                print ">>>>>>>>  Using CRUIJFF for systematics >>>>>>>>"
                systfunc = ROOT.TF1("systfunc",ROOT.cruijff, h.GetMean() - 2.3*h.GetRMS(),  h.GetMean() + 1.5*h.GetRMS(),5)
                systfunc.SetParameters(funct.GetParameter(0), funct.GetParameter(1), funct.GetParameter(2), 0., 0.) #15, 0.001)             
                systfunc.SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR")        
            systfunc.SetLineColor(ROOT.kRed)
            h.Fit("systfunc","M0R+")

        for par in range(funct.GetNpar()-1):
            pars.append(funct.GetParameter(par+1))
            if syst and ("Mean" in funct.GetParName(par+1) or "Sigma" in funct.GetParName(par+1)): 
                sys =1-systfunc.GetParameter(par+1)/funct.GetParameter(par+1)
                sys = sys*funct.GetParameter(par+1)
            else:
                sys = 0.
            errs.append(math.sqrt(sys*sys+funct.GetParError(par+1)*funct.GetParError(par+1)))
        chi2.append(funct.GetChisquare()/funct.GetNDF())

        h.SetTitle("Mass resolution for %d < m_{ll} <%d" %(mrange[i],mrange[i+1]))
        h.GetXaxis().SetTitle("m_{ll}^{RECO} / m_{ll}^{GEN} - 1")
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8)
        h.GetXaxis().SetRangeUser(fit_min,fit_max)

        h.Draw("E")
        funct.Draw("SAME")
#        systfunc.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextAlign(31)
        latex.SetTextSize(0.04)
        latex.SetNDC(True)
        latexCMS = ROOT.TLatex()
        latexCMS.SetTextFont(61)
        latexCMS.SetTextSize(0.055)
        latexCMS.SetNDC(True)
        latexCMSExtra = ROOT.TLatex()
        latexCMSExtra.SetTextFont(52)
        latexCMSExtra.SetTextSize(0.03)
        latexCMSExtra.SetNDC(True)

        latex.DrawLatex(0.95, 0.96, "(13 TeV)")

        cmsExtra = "Simulation" 
        latexCMS.DrawLatex(0.78,0.88,"CMS")
        yLabelPos = 0.84
        latexCMSExtra.DrawLatex(0.78,yLabelPos,"%s"%(cmsExtra))

        latexFit1 = ROOT.TLatex()
        latexFit1.SetTextFont(61)
        latexFit1.SetTextSize(0.035)
        latexFit1.SetNDC(True)
        latexFit1.DrawLatex(0.19, 0.84, "%d < m <%d" %(mrange[i],mrange[i+1]))
        
        latexFit = ROOT.TLatex()
        latexFit.SetTextFont(42)
        latexFit.SetTextSize(0.030)
        latexFit.SetNDC(True)        
        for par in range(funct.GetNpar()-1):
            yPos = 0.74-0.04*(float(par))
            latexFit.DrawLatex(0.19, yPos,"%s = %5.3g #pm %5.3g"%(funct.GetParName(par+1),funct.GetParameter(par+1),funct.GetParError(par+1)))
        latexFit.DrawLatex(0.19, 0.58, "#chi^{2}/ndf = %5.1f / %2.0f = %4.2f" %(funct.GetChisquare(),funct.GetNDF(),funct.GetChisquare()/funct.GetNDF()))
                
        saveas = "/MassRes_M%d_%d_%s_%s" %(mrange[i],mrange[i+1],rap,fit)
        c1.Print(output+saveas+".root")
        c1.Print(output+saveas+".C")
        c1.Print(output+saveas+".png")
        c1.Print(output+saveas+".pdf")
        
    print "DONE Fitting..."
    return pars,errs,chi2

def doFitWithSyst(hist,output,nrms,rapidity):
    print "######################################################"
    print "### FITTING HISTOS AND COMPUTING SYST  ERRORS      ###"
    print "######################################################"
    (sig     ,err,alp     ,aer,n     ,nerr) = doFit(hist,output,nrms,rapidity)
    (sig_down,_  ,alp_down,_  ,n_down,_   ) = doFit(hist,output,nrms*0.75,rapidity)
    (sig_up  ,_  ,alp_up  ,_  ,n_up  ,_   ) = doFit(hist,output,nrms*1.25,rapidity)
    
    for i in range(0,len(sig)):
        sys    = max(abs(1-sig_up[i]/sig[i]),abs(1-sig_down[i]/sig[i]))
        sys    = sys*sig[i]
        err[i] = math.sqrt(sys*sys+err[i]*err[i])

        sys    = max(abs(1-alp_up[i]/alp[i]),abs(1-alp_down[i]/alp[i]))
        sys    = sys*alp[i]
        aer[i] = math.sqrt(sys*sys+aer[i]*aer[i])

        sys     = max(abs(1-n_up[i]/n[i]),abs(1-n_up[i]/n[i]))
        sys     = sys*n[i]
        nerr[i] = math.sqrt(sys*sys+aer[i]*aer[i])

    print "############"
    print "### DONE ###"
    print "############"
    return sig,err,alp,aer,n,nerr
    

def drawMassResGeneric(hist,output,rapidity,funct="cruijff"):
    mass = []
    merr = []
    for i in range(len(mrange)-1):
        mass.append(mrange[i]+(mrange[i+1]-mrange[i])/2)
        merr.append((mrange[i+1]-mrange[i])/2)
    

    (pars,errs,chi2) = doFitGeneric(hist,output,rapidity,funct,True)
#    if "crystal" in funct: 
#        (pars2,_,_) = doFitGeneric(hist,output,rapidity,"cruijff")
#    else:
#        (pars2,_,_) = doFitGeneric(hist,output,rapidity,"crystal")
    
    c2 = ROOT.TCanvas("c2","c2",700,700)
    c2.cd()

    fun  = ROOT.TF1("fun","pol4")
    fun.SetParNames("A","B","C","D","E")            
    for i in range(fun.GetNpar()): 
        fun.ReleaseParameter(i)
        fun.SetParameter(i,0.)

    param = [ROOT.TGraphErrors(len(mass)) for x in range((len(pars)/len(mass))+1)] 
    res = ROOT.TGraphErrors(len(mass))
    for k,f in enumerate(param): 
        if k==4: 
            f.SetName("chi2")
        else : 
            f.SetName(hist[0].GetFunction(funct).GetParName(k+1))            

        for i in range(0,len(mass)):
            if k==4: 
                f.SetName("chi2")
                f.SetPoint(i,mass[i],chi2[i])
                f.SetPointError(i,merr[i],0)
            else: 
                f.SetPoint(i,mass[i],pars[i*4+k])
                f.SetPointError(i,merr[i],errs[i*4+k])
        
        if ("Sigma" in f.GetName()):
            res = param[k]

        f.SetMarkerStyle(20)
        f.SetMarkerSize(1.0)
        f.SetMarkerColor(ROOT.kBlue)
        f.SetLineColor(ROOT.kBlue)
        f.SetFillColor(0)
        f.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
        f.GetXaxis().SetRangeUser(mrange[0],mrange[len(mrange)-1])
        if ("chi2" in f.GetName()): 
            f.GetYaxis().SetRangeUser(0,20)            

        if "Sigma" in f.GetName():
            f.GetYaxis().SetRangeUser(0,0.15)

        if "AlphaR" in f.GetName():
            f.GetYaxis().SetRangeUser(0,.4)

        if "AlphaL" in f.GetName():
            f.GetYaxis().SetRangeUser(0.1,.6)

        if "Mean" in f.GetName():
            f.GetYaxis().SetRangeUser(-0.035,0.05)
                        
        f.Draw("AP E0")
        
        ## FIT PARAMETERS 
        for i in range(fun.GetNpar()): 
            fun.ReleaseParameter(i)
            fun.SetParameter(i,0.)
        
        if ("chi2" not in f.GetName()): 
            if ("Sigma" in f.GetName()):  
                print "Fitting Sigma"
                fun.SetParameters(0.,1E-5,-1.E-8,2E-12,-2E-16)
                fun.SetParLimits(1, 1.0E-6, 1.0E-4)
                fun.SetParLimits(2,-1.0E-7,-1.0E-9)
#                fun.SetParLimits(4,-3.0E-16,-1E-16)
                #                fun.FixParameter(3,0.)
#                fun.FixParameter(3,0.)
#                fun.FixParameter(4,0.)
            elif "AlphaR" in f.GetName(): 
                fun.SetParameters(0.25, 1E-6, -1.E-9, 1.E-12, -1.E-16)
                fun.SetParLimits(1, 1E-7, 1E-5)
                fun.SetParLimits(2, -1E-8 ,-1E-10)                
                fun.SetParLimits(3, 1E-14 ,1E-11)                
                fun.FixParameter(4,0.)
                fun.FixParameter(3,0.)
#                fun.FixParameter(2,0.)
#                fun.FixParameter(1,0.)
            elif "AlphaL" in f.GetName(): 
                fun.SetParameters(0.1,-1E-6, 1E-9, -1.E-13, 1E-16)
                fun.SetParLimits(1, -5E-5, -5E-7)
                fun.SetParLimits(2, 1E-10, 1E-8)                
                fun.SetParLimits(3, -1E-12, -1E-15)
#                fun.SetParLimits(4, 1E-18, 1E-10)
                fun.FixParameter(4,0.)
                fun.FixParameter(3,0.)
#                fun.FixParameter(2,0.)
            elif "Mean" in f.GetName():
                fun.SetParameters(0.004,-3E-5,1E-10,-1E-12,1.E-16)
                fun.SetParLimits(1,-1E-4,-1E-6)
                fun.SetParLimits(2, 1E-12,1E-8)
                fun.SetParLimits(3,-1E-12,-5E-14)
                fun.FixParameter(4,0.)
                
            f.Fit(fun,"MBFE+")            
            fun.Draw("SAME")

        
            latexFit = ROOT.TLatex()
            latexFit.SetTextFont(42)
            latexFit.SetTextSize(0.030)
            latexFit.SetNDC(True)        
            for par in range(fun.GetNpar()):
                yPos = 0.74-0.04*(float(par))
                latexFit.DrawLatex(0.19, yPos,"%s = %5.3g #pm %5.3g"%(fun.GetParName(par),fun.GetParameter(par),fun.GetParError(par)))
            latexFit.DrawLatex(0.19, 0.54, "#chi^{2}/ndf = %5.1f / %2.0f = %4.2f" %(fun.GetChisquare(),fun.GetNDF(),fun.GetChisquare()/fun.GetNDF()))
            
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextAlign(31)
        latex.SetTextSize(0.04)
        latex.SetNDC(True)
        latexCMS = ROOT.TLatex()
        latexCMS.SetTextFont(61)
        latexCMS.SetTextSize(0.055)
        latexCMS.SetNDC(True)
        latexCMSExtra = ROOT.TLatex()
        latexCMSExtra.SetTextFont(52)
        latexCMSExtra.SetTextSize(0.03)
        latexCMSExtra.SetNDC(True)
        
        latex.DrawLatex(0.95, 0.96, "(13 TeV)")
        
        cmsExtra = "Simulation" #splitline{Simulation}{Preliminary}"
        latexCMS.DrawLatex(0.19,0.88,"CMS")
        yLabelPos = 0.84
        latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))


        c2.SetGrid()
        saveas = "/%sVsMass_%s" %(f.GetName(),rapidity)
        c2.SaveAs(output+saveas+".png")
        c2.SaveAs(output+saveas+".pdf")
        c2.SaveAs(output+saveas+".root")
        c2.SaveAs(output+saveas+".C")
    
        ROOT.gPad.Update()
        c2.Clear()
    

    # PRINT FIT RESULTS!!!
    return res
    
    
def makeMassRes(inputfile,output,funct):
    style = setTDRStyle()
    ROOT.gStyle.SetTitleYOffset(1.45)
    ROOT.gStyle.SetTitleXOffset(1.45)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetStatX(.9)
    ROOT.gStyle.SetStatY(.9)
    
    hist_barrel = loadHistos(inputfile,"BB",1)
    hist_other  = loadHistos(inputfile,"BE",1)
    resBB  = drawMassResGeneric(hist_barrel,output,"BB",funct)
    resBE  = drawMassResGeneric(hist_other,output,"BE",funct)
    
    res = ROOT.TCanvas("res","res",700,700)
    res.cd()
    res.SetTickx()
    res.SetTicky()
    
    resBB.SetMarkerStyle(20)
    resBB.SetMarkerSize(1)
    resBB.SetMarkerColor(ROOT.kRed)
    resBB.SetLineColor(ROOT.kRed)
    resBB.SetFillColor(0)
    resBB.SetTitle("Dimuon mass resolution vs mass")
    resBB.GetYaxis().SetTitle("Dimuon Mass Resolution")
#    resBB.GetYaxis().SetTitleOffset(1.5)
    resBB.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    resBB.GetYaxis().SetRangeUser(0,.15)
    resBB.GetXaxis().SetRangeUser(mrange[0],mrange[len(mrange)-1])
    resBB.GetFunction("fun").SetLineColor(ROOT.kRed+1)
    resBB.Draw("AP E0")
    
    resBE.SetMarkerStyle(20)
    resBE.SetMarkerSize(1.0)
    resBE.SetMarkerColor(ROOT.kGreen+1)
    resBE.SetLineColor(ROOT.kGreen+1)
    resBE.SetFillColor(0)
    resBE.SetTitle("Dimuon mass resolution vs mass")
    resBE.GetYaxis().SetTitle("Dimuon Mass Resolution")
    resBE.GetYaxis().SetTitleOffset(1.5)
 #   resBE.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    resBE.GetYaxis().SetRangeUser(0,.15)
    resBE.GetXaxis().SetRangeUser(mrange[0],mrange[len(mrange)-1])
    resBE.GetFunction("fun").SetLineColor(ROOT.kGreen+2)
    resBE.Draw("PE0 SAME")

    latexFitBB = ROOT.TLatex()
    latexFitBB.SetTextFont(42)
    latexFitBB.SetTextSize(0.030)
    latexFitBB.SetNDC(True)        
    latexFitBB.SetTextColor(ROOT.kRed)

    latexFitBE = ROOT.TLatex()
    latexFitBE.SetTextFont(42)
    latexFitBE.SetTextSize(0.030)
    latexFitBE.SetNDC(True)        
    latexFitBE.SetTextColor(ROOT.kGreen+2)
    latexFitBB.DrawLatex(0.19, 0.78,"BB Category")
    latexFitBE.DrawLatex(0.60, 0.78,"BE+EE Category")
    for par in range(resBB.GetFunction("fun").GetNpar()):
        yPos = 0.74-0.04*(float(par))
        latexFitBB.DrawLatex(0.19, yPos,"%s = %5.3g #pm %5.3g"%(resBB.GetFunction("fun").GetParName(par),resBB.GetFunction("fun").GetParameter(par),resBB.GetFunction("fun").GetParError(par)))
        latexFitBE.DrawLatex(0.60, yPos,"%s = %5.3g #pm %5.3g"%(resBE.GetFunction("fun").GetParName(par),resBE.GetFunction("fun").GetParameter(par),resBE.GetFunction("fun").GetParError(par)))
    latexFitBB.DrawLatex(0.19, 0.54, "#chi^{2}/ndf = %5.1f / %2.0f = %4.2f" %(resBB.GetFunction("fun").GetChisquare(),resBB.GetFunction("fun").GetNDF(),resBB.GetFunction("fun").GetChisquare()/resBB.GetFunction("fun").GetNDF()))
    latexFitBE.DrawLatex(0.60, 0.54, "#chi^{2}/ndf = %5.1f / %2.0f = %4.2f" %(resBE.GetFunction("fun").GetChisquare(),resBE.GetFunction("fun").GetNDF(),resBE.GetFunction("fun").GetChisquare()/resBE.GetFunction("fun").GetNDF()))
        
#    leg = ROOT.TLegend(.35,.7,.50,.80,"","brNDC")
#    leg.AddEntry(resBB,"BB")
#    leg.AddEntry(resBE,"BE+EE")
#    leg.SetTextFont(42)
#    leg.SetBorderSize(0)
#    leg.SetTextSize(.04)
#    leg.Draw("SAME")

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.04)
    latex.SetNDC(True)
    latexCMS = ROOT.TLatex()
    latexCMS.SetTextFont(61)
    latexCMS.SetTextSize(0.055)
    latexCMS.SetNDC(True)
    latexCMSExtra = ROOT.TLatex()
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextSize(0.03)
    latexCMSExtra.SetNDC(True)
    
    latex.DrawLatex(0.95, 0.96, "(13 TeV)")
    
    cmsExtra = "Simulation" #splitline{Simulation}{Preliminary}"
    latexCMS.DrawLatex(0.19,0.88,"CMS")
    yLabelPos = 0.84
    latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
    
    res.SetGrid()
    saveas = "/MassResolutionVsMass" 
    res.SaveAs(output+saveas+".png")
    res.SaveAs(output+saveas+".pdf")
    res.SaveAs(output+saveas+".root")
    res.SaveAs(output+saveas+".C")
             
#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="makeMassRes.py [options]",description="Compute mass resolution",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i","--ifile", dest="inputfile",default="files/res_ZToMuMu_M_120_6000.root", help='Input filename')
    parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name to store results')
    parser.add_argument("-f","--funct",dest="funct", default="cruijff", help='function used')
    parser.add_argument("-x","--xrange", type=str, help='lower and upper x limit', nargs=1)
                    
#    parser.add_argument("-ncat","--ncategories", dest="ncat", type=int, default=3, help='number of categories')
    args = parser.parse_args()
    
    inputfile = args.inputfile
    output=args.output
#    ncat=args.ncat

    if args.xrange is not None:
        ranges = args.xrange[0].split(",")
        FITMIN = float(ranges[0])
        FITMAX = float(ranges[1])

    if not os.path.exists(args.output):
        os.makedirs(args.output);
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+args.output)
    
    print "Running on: %s " %(inputfile)
    print "Saving result in: %s" %(output)

    makeMassRes(inputfile,output,args.funct)
    print "DONE"
