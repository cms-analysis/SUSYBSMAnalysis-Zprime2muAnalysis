import FWCore.ParameterSet.Config as cms

# Hack to set the python path for now in 1_6_X.
import sys, os
sys.path.insert(0, os.path.join(os.environ['CMSSW_BASE'],
                                'src/SUSYBSMAnalysis/Zprime2muAnalysis/data'))

from Zprime2muAnalysisCommon_cff import *
from Zprime2muResolution_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_1.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_2.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_3.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_4.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_5.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_6.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_7.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_8.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_9.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_10.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_11.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_12.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_13.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_14.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_15.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_16.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_17.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_18.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_19.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_0_6/ZPSSMmumu_M1000_Mcut400/reco_20.root'
    ]

process = makeZprime2muAnalysisProcess(files) #, performTrackReReco=True)
#selectAlignment(process)
attachResolution(process)

process.Zprime2muResolution.verbosity = 2
process.Zprime2muResolution.dateHistograms = False

#print process.dumpConfig()

# Example of how to replace one of the parameter sets:
'''
process.Zprime2muResolution.exampleSet = cms.PSet(
    peakMass = cms.double(200),
    lowerMassWin = cms.double(0),
    upperMassWin = cms.double(300),
    binSize = cms.int32(25),
    maxTrigMass = cms.double(0.3))
process.Zprime2muResolution.dataSet = 'exampleSet'
'''
