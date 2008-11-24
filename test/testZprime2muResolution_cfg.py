import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muResolution_cfi import *

# Use a RelVal Z0->mumu sample for now until official samples are available.
files = [
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/12492C75-CE78-DD11-8363-001617C3B5F4.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/1C40E3E6-CC78-DD11-980F-001617C3B64C.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/1C66C875-CE78-DD11-B6D7-000423D6CA6E.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/221D6ED5-CF78-DD11-BC9D-000423D996C8.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/548D5E2C-D778-DD11-A8F1-001D09F24259.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/5E3DFEDA-CE78-DD11-915F-000423D98DB4.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/80D656E6-D778-DD11-973E-001617E30CA4.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/868E8D72-CF78-DD11-9418-001617C3B6CE.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/88F58F12-CE78-DD11-BAB9-000423D985E4.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/A6EE09FF-CD78-DD11-A2CB-000423D98804.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/AE99746B-D078-DD11-ADD6-001617DBD230.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/BAB35A53-D778-DD11-A88A-001D09F28F11.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/D204572A-D778-DD11-951E-0030487C5CFA.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/D4E91504-D078-DD11-B961-000423D9870C.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/D8FC112C-CC78-DD11-9DD4-0019DB2F3F9B.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0000/DED1077D-CD78-DD11-AC47-001617DF785A.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/AAA5B85B-F878-DD11-8CDF-000423D6BA18.root',
    '/store/relval/CMSSW_2_1_6/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/B417007E-EA78-DD11-86CC-001617C3B77C.root'
]

#files = poolAllFiles('/scratchdisk3/tucker/CMSSW_2_1_6/RelValZMM/*.root')

process = makeZprime2muAnalysisProcess(files, useTrigger=True)
attachResolution(process)

process.Zprime2muResolution.verbosity = 2
process.Zprime2muResolution.dateHistograms = False

# Example of how to replace one of the parameter sets (e.g. for
# running on the Z0->mumu sample above):
process.Zprime2muResolution.Z0params = cms.PSet(
    peakMass = cms.double(91.2),
    lowerMassWin = cms.double(0),
    upperMassWin = cms.double(300),
    binSize = cms.int32(25),
    maxTrigMass = cms.double(0.3)
    )
process.Zprime2muResolution.dataSet = 'Z0params'

#print process.dumpConfig()

