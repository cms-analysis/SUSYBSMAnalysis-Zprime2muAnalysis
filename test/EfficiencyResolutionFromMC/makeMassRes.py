#!/usr/bin/python

# import ROOT in batch mode
import sys
import argparse
import math

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv
    
def loadHistos(inputfile,region,rebin):
    _file = ROOT.TFile(inputfile)

    histos = []

    Res_0to500     = ROOT.TH1F()
    Res_500to1000  = ROOT.TH1F()
    Res_1000to1500 = ROOT.TH1F()
    Res_1500to2000 = ROOT.TH1F()
    Res_2000to2500 = ROOT.TH1F()
    Res_2500to3000 = ROOT.TH1F()
    Res_3000to3500 = ROOT.TH1F()
    Res_3500to4000 = ROOT.TH1F()
    Res_4000to4500 = ROOT.TH1F()
    Res_4500to5000 = ROOT.TH1F()
    
    Res_0to500     .SetDirectory(0)
    Res_500to1000  .SetDirectory(0)
    Res_1000to1500 .SetDirectory(0)
    Res_1500to2000 .SetDirectory(0)
    Res_2000to2500 .SetDirectory(0)
    Res_2500to3000 .SetDirectory(0)
    Res_3000to3500 .SetDirectory(0)
    Res_3500to4000 .SetDirectory(0)
    Res_4000to4500 .SetDirectory(0)
    Res_4500to5000 .SetDirectory(0)
    ROOT.TH1.AddDirectory(ROOT.kFALSE)
    
    if ("barrel" in region):
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500BB").Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000BB").Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500BB").Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000BB").Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500BB").Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000BB").Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500BB").Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000BB").Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500BB").Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000BB").Clone()
    elif ("other" in region):
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500BEp").Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000BEp").Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500BEp").Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000BEp").Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500BEp").Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000BEp").Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500BEp").Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000BEp").Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500BEp").Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000BEp").Clone()

        Res_0to500     .Add( _file.Get("DileptonMassResVMass_0to500BEm").Clone()    )
        Res_500to1000  .Add( _file.Get("DileptonMassResVMass_500to1000BEm").Clone() )
        Res_1000to1500 .Add( _file.Get("DileptonMassResVMass_1000to1500BEm").Clone())
        Res_1500to2000 .Add( _file.Get("DileptonMassResVMass_1500to2000BEm").Clone())
        Res_2000to2500 .Add( _file.Get("DileptonMassResVMass_2000to2500BEm").Clone())
        Res_2500to3000 .Add( _file.Get("DileptonMassResVMass_2500to3000BEm").Clone())
        Res_3000to3500 .Add( _file.Get("DileptonMassResVMass_3000to3500BEm").Clone())
        Res_3500to4000 .Add( _file.Get("DileptonMassResVMass_3500to4000BEm").Clone())
        Res_4000to4500 .Add( _file.Get("DileptonMassResVMass_4000to4500BEm").Clone())
        Res_4500to5000 .Add( _file.Get("DileptonMassResVMass_4500to5000BEm").Clone())
                    
    elif ("all" in region):
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500").Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000").Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500").Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000").Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500").Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000").Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500").Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000").Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500").Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000").Clone()
        
    elif ("BB" in region):
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500BB").Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000BB").Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500BB").Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000BB").Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500BB").Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000BB").Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500BB").Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000BB").Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500BB").Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000BB").Clone()
    elif ("BEp" in region):
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500BEp").Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000BEp").Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500BEp").Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000BEp").Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500BEp").Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000BEp").Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500BEp").Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000BEp").Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500BEp").Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000BEp").Clone()
    elif ("BEm" in region):
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500BEm").Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000BEm").Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500BEm").Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000BEm").Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500BEm").Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000BEm").Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500BEm").Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000BEm").Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500BEm").Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000BEm").Clone()
    else:
        Res_0to500     = _file.Get("DileptonMassResVMass_0to500%s"%region).Clone()
        Res_500to1000  = _file.Get("DileptonMassResVMass_500to1000%s"%region).Clone()
        Res_1000to1500 = _file.Get("DileptonMassResVMass_1000to1500%s"%region).Clone()
        Res_1500to2000 = _file.Get("DileptonMassResVMass_1500to2000%s"%region).Clone()
        Res_2000to2500 = _file.Get("DileptonMassResVMass_2000to2500%s"%region).Clone()
        Res_2500to3000 = _file.Get("DileptonMassResVMass_2500to3000%s"%region).Clone()
        Res_3000to3500 = _file.Get("DileptonMassResVMass_3000to3500%s"%region).Clone()
        Res_3500to4000 = _file.Get("DileptonMassResVMass_3500to4000%s"%region).Clone()
        Res_4000to4500 = _file.Get("DileptonMassResVMass_4000to4500%s"%region).Clone()
        Res_4500to5000 = _file.Get("DileptonMassResVMass_4500to5000%s"%region).Clone()
        
    histos.append(Res_0to500    )
    histos.append(Res_500to1000 )
    histos.append(Res_1000to1500)
    histos.append(Res_1500to2000)
    histos.append(Res_2000to2500)
    histos.append(Res_2500to3000)
    histos.append(Res_3000to3500)
    histos.append(Res_3500to4000)
    histos.append(Res_4000to4500)
    histos.append(Res_4500to5000)

#        if(h.Integral() < 5000.): 
#            print "Rebinning histos with value: %d" %(rebin*2)
#            h.Rebin(rebin*2) 
#        else:
#            print "Rebinning histos with value: %d" %(rebin)            
#            h.Rebin(rebin) 
        
    _file.Close()
    return histos

def doFit(hist,output,nrms,rap="BB"):
    c1 = ROOT.TCanvas("c1","c1",700,700)
    c1.cd()

    leg = ROOT.TLegend(.20,.7,.30,.80,"","brNDC");
    leg.SetTextFont(42);
    leg.SetBorderSize(0);
    leg.SetTextSize(.02);

    mrange = [0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]
    sig    = []
    sige   = []
    alpha  = []
    alphae = []
    n      = []
    ne     = []
    for i,h in enumerate(hist):
        print "+++++++++++++++++++++++++++++++++++++++++"
        print "Fitting histogram for %d < m_{ll} <%d" %(mrange[i],mrange[i+1])
        print "+++++++++++++++++++++++++++++++++++++++++\n"

        fit_min = h.GetMean() - 1.3*nrms*h.GetRMS() 
        fit_max = h.GetMean() + 0.7*nrms*h.GetRMS()

        # fit with a gaussian to use parameters of the fit for the CB...
        gaus = ROOT.TF1("gaus_%s"%(nrms),"gaus",fit_min,fit_max)
        gaus.SetLineColor(ROOT.kBlue)
        gaus.SetParameters(0,h.GetMean(),h.GetRMS())
        h.Fit("gaus_%s"%(nrms),"M0R+")
        
        crystal = ROOT.TF1("crystal_%s"%(nrms),"crystalball",fit_min,fit_max)
        crystal.SetLineColor(ROOT.kRed)
        
        tmp_mean = gaus.GetParameter(1)
        tmp_sig  = gaus.GetParameter(2)
        min_mean = tmp_mean+0.5*tmp_mean
        max_mean = tmp_mean-0.5*tmp_mean
        crystal.SetParameters(gaus.GetParameter(0), tmp_mean, tmp_sig, 1.4, 2)
#        crystal.SetParLimits(0, 0, 2*gaus.GetParameter(0)) # //const
        crystal.SetParLimits(1, min_mean, max_mean)
        crystal.SetParLimits(2, 0, 2.5*h.GetRMS())# //sigma
        crystal.SetParLimits(3, 0.5, 2.)# //alpha
#        crystal.SetParLimits(4, 0., 3.)# //alpha
        
#        crystal.SetParLimits(4, 0, 5.)
#            crystal.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1), gaus.GetParameter(2), -1.6, 7.85)
#            crystal.SetParLimits(3, -3, 3)# //alpha
        
        h.Fit("crystal_%s"%(nrms),"M0R+")
        
        # crystal
        sig   .append(crystal.GetParameter(2))
        sige  .append(crystal.GetParError(2))
        alpha .append(crystal.GetParameter(3))  
        alphae.append(crystal.GetParError(3))
        n     .append(crystal.GetParameter(4))  #exponential parameter 
        ne    .append(crystal.GetParError(4))   #exponential parameter

        h.SetTitle("Mass resolution for %d < m_{ll} <%d" %(mrange[i],mrange[i+1]))
        h.GetXaxis().SetTitle("m_{ll}^{RECO} / m_{ll}^{GEN} - 1")
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.7)
    
        if (i==0):
            leg.AddEntry(h,"DY simulation")
#            leg.AddEntry(h.GetFunction("gaus_%s"%(nrms)),"Gaussian Fit","L")
            leg.AddEntry(h.GetFunction("crystal_%s"%(nrms)),"CrystalBall Fit","L")
            
        h.Draw("E")
        crystal.Draw("SAME")
        leg.Draw("SAME")
        
        ROOT.gSystem.MakeDirectory("%s/%1.1fRMS" %(output,nrms))
        saveas = "/%1.1fRMS/MassRes_M%d_%d_%s" %(nrms,mrange[i],mrange[i+1],rap)
        c1.SaveAs(output+saveas+".root")
        c1.SaveAs(output+saveas+".C")
        c1.SaveAs(output+saveas+".png")
        c1.SaveAs(output+saveas+".pdf")

    print "DONE Fitting..."
    return sig,sige,alpha,alphae,n,ne

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

    
def drawMassRes(hist,output,rapidity,nrms):
    mass = [250,750,1250,1750,2250,2750,3250,3750,4250,4750]
    merr = [250,250,250,250,250,250,250,250,250,250]

#    (c_sig,c_err,g_sig,g_err) = doFit(hist,output,nrms,rapidity,2)
    (sig,err,alpha,alphae,n,nerr) = doFitWithSyst(hist,output,nrms,rapidity)
    
    c2 = ROOT.TCanvas("c2","c2",700,700)
    c2.cd()

    fun = ROOT.TF1("fun","pol2")
    fun.SetParameters(0.,0.,0.)
  
    res_crystal  = ROOT.TGraphErrors(10)
    res_crystal.SetName("crystal")
    for i in range(0,len(mass)):
        res_crystal.SetPoint(i,mass[i],sig[i])
        res_crystal.SetPointError(i,merr[i],err[i])

    res_crystal.SetMarkerStyle(22)
    res_crystal.SetMarkerColor(ROOT.kRed)
    res_crystal.SetLineColor(ROOT.kRed)
    res_crystal.SetFillColor(0)
    res_crystal.SetTitle("Dimuon mass resolution vs mass")
    res_crystal.GetYaxis().SetTitle("Dimuon Mass Resolution")
    res_crystal.GetYaxis().SetTitleOffset(1.5)
    res_crystal.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    res_crystal.GetYaxis().SetRangeUser(0,.1)
    res_crystal.GetXaxis().SetRangeUser(0,5000)
    fun.SetParameters(0.,0.,0.)
    res_crystal.Fit(fun,"M+")
    res_crystal.GetFunction("fun").SetLineColor(ROOT.kRed+2)
    res_crystal.Draw("AP E0")

    leg = ROOT.TLegend(.35,.7,.50,.80,"","brNDC")
    leg.AddEntry(res_crystal,"CrystalBall Fit")
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.SetTextSize(.02)
    leg.Draw("SAME")

    c2.SetGrid()
    saveas = "/%1.1fRMS/MassResolutionVsMass_%s" %(nrms,rapidity)
    c2.SaveAs(output+saveas+".png")
    c2.SaveAs(output+saveas+".pdf")
    c2.SaveAs(output+saveas+".root")
    c2.SaveAs(output+saveas+".C")
    
    ROOT.gPad.Update()
    c2.Clear()
    
    alp_crystal  = ROOT.TGraphErrors(10)
    alp_crystal.SetName("crystal")
    for i in range(0,len(mass)):
        alp_crystal.SetPoint(i,mass[i],alpha[i])
        alp_crystal.SetPointError(i,merr[i],alphae[i])

    alp_crystal.SetMarkerStyle(22)
    alp_crystal.SetMarkerColor(ROOT.kRed)
    alp_crystal.SetLineColor(ROOT.kRed)
    alp_crystal.SetFillColor(0)
    alp_crystal.GetYaxis().SetRangeUser(0,2.)
    alp_crystal.GetXaxis().SetRangeUser(0,5000)
    alp_crystal.SetTitle("CB parameter alpha vs mass")
    alp_crystal.GetYaxis().SetTitle("CB parameter alpha")
    alp_crystal.GetYaxis().SetTitleOffset(1.5)
    alp_crystal.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    fun.SetParameters(0.,0.,0.)
    fun.FixParameter(1,0.)
    fun.FixParameter(2,0.)
    alp_crystal.Fit(fun,"M+")
    alp_crystal.GetFunction("fun").SetLineColor(ROOT.kRed+2)
    alp_crystal.Draw("AP E0")

    saveas = "/%1.1fRMS/AlphaVsMass_%s" %(nrms,rapidity)
    c2.SaveAs(output+saveas+".png")
    c2.SaveAs(output+saveas+".pdf")
    c2.SaveAs(output+saveas+".root")
    
    n_crystal  = ROOT.TGraphErrors(10)
    n_crystal.SetName("crystal")
    for i in range(0,len(mass)):
        n_crystal.SetPoint(i,mass[i],n[i])
        n_crystal.SetPointError(i,merr[i],nerr[i])

    n_crystal.SetMarkerStyle(22)
    n_crystal.SetMarkerColor(ROOT.kRed)
    n_crystal.SetLineColor(ROOT.kRed)
    n_crystal.SetFillColor(0)
    n_crystal.GetYaxis().SetRangeUser(0,3.)
    n_crystal.GetXaxis().SetRangeUser(0,5000)
    n_crystal.SetTitle("Parametrization of n vs mass")
    n_crystal.GetYaxis().SetTitle("CB parameter n")
    n_crystal.GetYaxis().SetTitleOffset(1.5)
    n_crystal.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    fun.SetParameters(0.,0.,0.)
    fun.FixParameter(1,0.)
    fun.FixParameter(2,0.)
    n_crystal.Fit(fun,"M+")
    n_crystal.GetFunction("fun").SetLineColor(ROOT.kRed+2)
    n_crystal.Draw("AP E0")

    saveas = "/%1.1fRMS/BremsVsMass_%s" %(nrms,rapidity)
    c2.SaveAs(output+saveas+".png")
    c2.SaveAs(output+saveas+".pdf")
    c2.SaveAs(output+saveas+".root")
    
    
    # PRINT FIT RESULTS!!!
    mrange = [0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]
    ndf = []
    chi = []
    for h in hist:
        chi.append(h.GetFunction("crystal_%s"%nrms).GetChisquare())
        ndf.append(h.GetFunction("crystal_%s"%nrms).GetNDF())

    print "|---------------------------------------------------------------------------------------|"
    print "|                              CRYSTAL BALL PARAMETRIZATION                             |"
    print "|---------------------------------------------------------------------------------------|"
    print "|      mll  %s     | Chi2/n.d.f. |   Sigma  [%%]   |      Alpha      |       N          |" %(rapidity)
    print "|---------------------------------------------------------------------------------------|" 
    for i in range(0,len(mass)):
        print "| %4d < mll < %4d | %6.1f / %2.0f | %5.3f +/- %5.3f | %5.3f +/- %5.3f | %5.3f +/- %5.3f |" %(mass[i]-merr[i], mass[i]+merr[i], chi[i], ndf[i],
	        sig[i]*100, err[i]*100, alpha[i], alphae[i], n[i], nerr[i])
    print "|---------------------------------------------------------------------------------------------------------|" 
    print "| res(m)   | (%.2e +/- %.1e) + (%.2e +/- %.1e) x + (% .2e +/- %.1e) x^2 | %5.2f / %2.0f |" %(res_crystal.GetFunction("fun").GetParameter(0), res_crystal.GetFunction("fun").GetParError(0),
          res_crystal.GetFunction("fun").GetParameter(1), res_crystal.GetFunction("fun").GetParError(1),
          res_crystal.GetFunction("fun").GetParameter(2), res_crystal.GetFunction("fun").GetParError(2),
          res_crystal.GetFunction("fun").GetChisquare(),res_crystal.GetFunction("fun").GetNDF())
    print "| alpha(m) | (%.2e +/- %.1e) + (%.2e +/- %.1e) x + (% .2e +/- %.1e) x^2 | %5.2f / %2.0f |" %(alp_crystal.GetFunction("fun").GetParameter(0), alp_crystal.GetFunction("fun").GetParError(0),
          alp_crystal.GetFunction("fun").GetParameter(1), alp_crystal.GetFunction("fun").GetParError(1),
          alp_crystal.GetFunction("fun").GetParameter(2), alp_crystal.GetFunction("fun").GetParError(2),
          alp_crystal.GetFunction("fun").GetChisquare(),alp_crystal.GetFunction("fun").GetNDF())
    print "| n(m)     | (%.2e +/- %.1e) + (%.2e +/- %.1e) x + (% .2e +/- %.1e) x^2 | %5.2f / %2.0f |"  %(n_crystal.GetFunction("fun").GetParameter(0), n_crystal.GetFunction("fun").GetParError(0),
          n_crystal.GetFunction("fun").GetParameter(1), n_crystal.GetFunction("fun").GetParError(1),
          n_crystal.GetFunction("fun").GetParameter(2), n_crystal.GetFunction("fun").GetParError(2),
          n_crystal.GetFunction("fun").GetChisquare(),n_crystal.GetFunction("fun").GetNDF())      
    print "|---------------------------------------------------------------------------------------------------------|" 
    
    return res_crystal
    
    
def makeMassRes(inputfile,output,nrms,ncat):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetStatX(.9)
    ROOT.gStyle.SetStatY(.9)

    ROOT.gSystem.MakeDirectory(output)
    
    
    if (ncat==2): 
        hist_barrel = loadHistos(inputfile,"barrel",1)
        hist_other  = loadHistos(inputfile,"other",1)
        hist_all    = loadHistos(inputfile,"all",1)
        drawMassRes(hist_barrel,output,"barrel",nrms)
        drawMassRes(hist_other,output,"other",nrms)
        drawMassRes(hist_all,output,"all",nrms)    
    elif (ncat==3):
         hist_BB  = loadHistos(inputfile,"BB",1)
         hist_BEp = loadHistos(inputfile,"BEp",1)
         hist_BEm = loadHistos(inputfile,"BEm",1)
         hist_all = loadHistos(inputfile,"all",1)
         resBB  = drawMassRes(hist_BB,  output, "BB" , nrms)
         resBEp = drawMassRes(hist_BEp, output, "BEp", nrms)
         resBEm = drawMassRes(hist_BEm, output, "BEm", nrms)
         resAll = drawMassRes(hist_all, output, "all", nrms)    
    
    print resBB
    print resBEp
    print resBEm
    
    res = ROOT.TCanvas("res","res",700,700)
    res.cd()
    res.SetTickx()
    res.SetTicky()
    
    resBB.SetMarkerStyle(22)
    resBB.SetMarkerColor(ROOT.kRed)
    resBB.SetLineColor(ROOT.kRed)
    resBB.SetFillColor(0)
    resBB.SetTitle("Dimuon mass resolution vs mass")
    resBB.GetYaxis().SetTitle("Dimuon Mass Resolution")
    resBB.GetYaxis().SetTitleOffset(1.5)
    resBB.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    resBB.GetYaxis().SetRangeUser(0,.2)
    resBB.GetXaxis().SetRangeUser(0,5000)
    resBB.GetFunction("fun").SetLineColor(ROOT.kRed+1)
    resBB.Draw("AP E0")
    
         
    resBEp.SetMarkerStyle(22)
    resBEp.SetMarkerColor(ROOT.kGreen+1)
    resBEp.SetLineColor(ROOT.kGreen+1)
    resBEp.SetFillColor(0)
    resBEp.SetTitle("Dimuon mass resolution vs mass")
    resBEp.GetYaxis().SetTitle("Dimuon Mass Resolution")
    resBEp.GetYaxis().SetTitleOffset(1.5)
    resBEp.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    resBEp.GetYaxis().SetRangeUser(0,.2)
    resBEp.GetXaxis().SetRangeUser(0,5000)
    resBEp.GetFunction("fun").SetLineColor(ROOT.kGreen+2)
    resBEp.Draw("PE0 SAME")
    
    resBEm.SetMarkerStyle(22)
    resBEm.SetMarkerColor(ROOT.kBlue+1)
    resBEm.SetLineColor(ROOT.kBlue+1)
    resBEm.SetFillColor(0)
    resBEm.SetTitle("Dimuon mass resolution vs mass")
    resBEm.GetYaxis().SetTitle("Dimuon Mass Resolution")
    resBEm.GetYaxis().SetTitleOffset(1.5)
    resBEm.GetXaxis().SetTitle("m(#mu^{+}#mu^{-}) [GeV]")
    resBEm.GetYaxis().SetRangeUser(0,.2)
    resBEm.GetXaxis().SetRangeUser(0,5000)
    resBEm.GetFunction("fun").SetLineColor(ROOT.kBlue+2)
    resBEm.Draw("PE0 SAME")
        
    leg = ROOT.TLegend(.35,.7,.50,.80,"","brNDC")
    leg.AddEntry(resBB,"BB")
    leg.AddEntry(resBEp,"BEp")
    leg.AddEntry(resBEm,"BEm")                
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.SetTextSize(.02)
    leg.Draw("SAME")
    
    res.SetGrid()
    saveas = "/%1.1fRMS/MassResolutionVsMass_3CAT" %(nrms)
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
    parser.add_argument("-n","--nrms",dest="nrms", type=float, default=2.0, help='number of RMS used in the fit')
    parser.add_argument("-ncat","--ncategories", dest="ncat", type=int, default=3, help='number of categories')
    args = parser.parse_args()
    
    inputfile = args.inputfile
    nrms= args.nrms
    output=args.output
    ncat=args.ncat
    
    print "Running on: %s with %d RMS and %d categories" %(inputfile,nrms,ncat)
    print "Saving result in: %s" %(output)

    makeMassRes(inputfile,output,nrms, ncat)
    print "DONE"
