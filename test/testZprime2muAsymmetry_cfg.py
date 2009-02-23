import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAsymmetry_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/ZPSSMmumu_M1000_Mcut400_10TeV_IDEAL_V9_RAW2DIGI_RECO_1-10.root'
]

process = makeZprime2muAnalysisProcess(files)
attachAsymmetry(process)

#process.Zprime2muAsymmetry.numFits = 4
#process.Zprime2muAsymmetry.useCachedParams = True
#process.Zprime2muAsymmetry.dataSet = 'Z0'
#process.Zprime2muAsymmetry.Z0 = cms.PSet(
#    genSampleFiles = cms.vstring('/scratchdisk3/tucker/CMSSW_2_1_6/RelValZMM/1C66C875-CE78-DD11-B6D7-000423D6CA6E.root'),
#    peakMass = cms.double(91.2), # no peak, but set this for completeness
#    massDistType = cms.int32(3),
#    maxPt = cms.double(200.0),
#    maxRapidity = cms.double(3.5),
#    outputFileBase = cms.untracked.string('Z0'),
#    genWindow = cms.vdouble(80.0, 100.0),
#    fitWindowOnPeak = cms.vdouble(80.0, 100.0),
#    fitWindowOffPeak = cms.vdouble(200.0, 500.0),
#    #   order: cos_cs, rapidity, pT, phi, mass, phi_cs
#    recSigma = cms.vdouble(5.8e-4, 8.6e-3, 3.9, 8.6e-2, 6.14, 9.5e-2)
#)
