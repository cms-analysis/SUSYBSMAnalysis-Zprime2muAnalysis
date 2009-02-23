import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_2_1_12/ZPSSMmumu_M1000_Mcut400_10TeV_IDEAL_V9_RAW2DIGI_RECO_1-10.root'
    ]

process = makeZprime2muAnalysisProcess(files)
attachAnalysis(process, 'Zprime2muExample',
               verbosity = 2,
               nBins     = cms.int32(40),
               lowerMass = cms.double(400),
               upperMass = cms.double(2000)
               )
