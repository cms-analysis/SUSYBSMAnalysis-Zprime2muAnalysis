#!/usr/bin/env python
'''
- Running module as stand alone code
'''
import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import *
import bckgshape_syst
import bckgshape_tools
#hists_dir = '/afs/cern.ch/work/c/cschnaib/Zprime2muAnalysis/DataMCSpectraComparision/mc/76X_v2/' # should be here oops
hists_dir = '/afs/cern.ch/work/c/cschnaib/DataMCSpectraComparision/mc/76X_v2/'

set_zp2mu_style()
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.TH1.AddDirectory(0)
#ps = plot_saver('plots/shapesyst')

rebin = 40
low = 40
fitlow = 400
#fithigh = 5500
fithigh = 5500
high = 5500
# Run2015(B+C+D)
int_lumi = 2800.

# Make nominal MC histogram
histMC = bckgshape_tools.make_mc_hist(True,False,int_lumi,low,high,rebin,hists_dir)
# Make MC histogram with PI background
histMCPI = bckgshape_tools.make_mc_hist(True,True,int_lumi,low,high,rebin,hists_dir)
# Make nominal MC fit and saves plot
nominal = bckgshape_tools.fit_hist(histMC,'nominal',low,high,fitlow,fithigh,rebin)
nominalPI = bckgshape_tools.fit_hist_PI(histMCPI,'nominalPI',low,high,fitlow,fithigh,rebin)
# Plot Bckg and Bckg+PI fit functions
bckgshape_tools.plot_two(nominal,'Nominal',histMC,nominalPI,'Nominal+PI',histMCPI,'nominal_PI','func')
# 15%
bckgshape_tools.fit_syst_func(nominal,histMC,bckgshape_syst.unc15,'unc15',rebin,low,high,fitlow,fithigh,'15')
# 20%
bckgshape_tools.fit_syst_func(nominal,histMC,bckgshape_syst.unc20,'unc20',rebin,low,high,fitlow,fithigh,'20')
# 25%
bckgshape_tools.fit_syst_func(nominal,histMC,bckgshape_syst.unc25,'unc25',rebin,low,high,fitlow,fithigh,'25')
# 50%
bckgshape_tools.fit_syst_func(nominal,histMC,bckgshape_syst.unc50,'unc50',rebin,low,high,fitlow,fithigh,'50')
print 'end'
