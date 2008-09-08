import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muMassReach_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/dy_above400_RECO.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/dy_above200_RECO.root'
#   'file:/data0/slava/data/ZPSSMmumu_M1000_Mcut400_1.root',
#   'file:/data0/slava/data/ZPSSMmumu_M1000_Mcut400_1.root',
#   'file:/data0/slava/data/ZPSSMmumu_M1000_Mcut400_1.root'
    ]

process = makeZprime2muAnalysisProcess(files)

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
#process.Zprime2muMassReach.intLumi       = 0.34

