import FWCore.ParameterSet.Config as cms

# Hack to set the python path for now in 1_6_X.
import sys, os
sys.path.insert(0, os.path.join(os.environ['CMSSW_BASE'],
                                'src/SUSYBSMAnalysis/Zprime2muAnalysis/data'))

from Zprime2muAnalysisCommon_cff import *
from Zprime2muResolution_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root'
#    'file:/scratchdisk2/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root'
    ]

process = makeZprime2muAnalysisProcess(files)
#   , doingElectrons=True,
#   flavorsForDileptons=diElectrons)

attachResolution(process)

#process.Zprime2muResolution.verbosity = 2
#process.Zprime2muResolution.dateHistograms = False

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
