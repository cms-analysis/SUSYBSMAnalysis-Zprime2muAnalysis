import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muResolution_cfi import *

files = [
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/0E637D36-9999-DD11-B684-001617E30D52.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/1E1BE784-A099-DD11-BCD6-000423D94E1C.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/26073DBE-9F99-DD11-B5E8-001617E30D40.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/2EA99A9A-A299-DD11-A6A3-000423D98930.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/3C4DCDF8-9699-DD11-9812-0016177CA778.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/4403B8BA-9F99-DD11-819D-001617C3B6E8.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/4EC45CE9-9799-DD11-B701-001617C3B5D8.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/5A8FBA47-A099-DD11-9920-000423D9A212.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/7437A584-9899-DD11-8833-001617E30D4A.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/76A00883-9999-DD11-AA86-000423D94AA8.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/80EBC746-FD99-DD11-AC24-000423D952C0.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/98F25C73-9F99-DD11-BBE8-001617E30D00.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/E8F9D8E9-A099-DD11-9F96-000423D98FBC.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/F2BB9497-9D99-DD11-A514-001617E30D12.root',
    '/store/relval/CMSSW_2_1_10/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/FA36E8F6-9E99-DD11-9161-000423D98B5C.root'
    ]

process = makeZprime2muAnalysisProcess(files)
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

