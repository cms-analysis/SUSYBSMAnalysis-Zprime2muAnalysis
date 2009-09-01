import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muMassReach_cfi import *

gtag = "MC_31X_V3"
dir = "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/"
attr = "_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_"

signame = "ZPSSMmumu_M1000_Mcut400"
sigdir = signame + "-" + gtag + "/"
signal = "PYTHIA6_" + signame + attr + gtag

bckg1name = "DYmumu_Mcut500"
bckg2name = "DYmumu_Mcut200"

bckg1dir = bckg1name + "-" + gtag + '/'
bckg2dir = bckg2name + "-" + gtag + '/'

bckg1 = 'PYTHIA6_' + bckg1name + attr + gtag
bckg2 = 'PYTHIA6_' + bckg2name + attr + gtag

files = [
    dir + sigdir + signal + "_1.root",
    dir + bckg1dir + bckg1 + "_1.root",
    dir + bckg2dir + bckg2 + "_1.root"
#    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root',
#    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/dy_above400_RECO.root',
#    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/dy_above200_RECO.root'
    ]

# process = makeZprime2muAnalysisProcess(files)
process = makeZprime2muAnalysisProcess(files, maxEvents=-1, skipPAT=True, disableElectrons=True, conditionsGlobalTag='STARTUP31X_V2::All') #, useTrigger=False)

#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/dy_for1TeV_official.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zssm1000_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi1000_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zssm1500_private_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi1500_private_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi1500_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zssm2000_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi2000_massreach.cfi")
    
attachMassReach(process)
process.Zprime2muMassReach.verbosity = 0

#process.Zprime2muMassReach.dataSet       = "Zssm1200_EPE"

#process.Zprime2muMassReach.dataSet       = "Zpsi1000"
#process.Zprime2muMassReach.dataSet       = "Zssm1500"
#process.Zprime2muMassReach.dataSet       = "Zpsi1500"
#process.Zprime2muMassReach.dataSet       = "Zssm2000"
#process.Zprime2muMassReach.dataSet       = "Zpsi2000"
#process.Zprime2muMassReach.DYEvents      = True
#process.Zprime2muMassReach.BackgroundFit = True
#process.Zprime2muMassReach.FixedMass     = False
#process.Zprime2muMassReach.FixedFWHM     = False
#process.Zprime2muMassReach.ExpPlots      = True
#process.Zprime2muMassReach.BinnedFit     = True
#process.Zprime2muMassReach.intLumi       = 0.05
#process.Zprime2muMassReach.intLumi       = 0.1

