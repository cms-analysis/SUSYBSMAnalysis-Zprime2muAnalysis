import FWCore.ParameterSet.Config as cms

# Hack to set the python path for now in 1_6_X.
import sys, os
sys.path.insert(0, os.path.join(os.environ['CMSSW_BASE'],
                                'src/SUSYBSMAnalysis/Zprime2muAnalysis/data'))

from Zprime2muAnalysisCommon_cff import *

files = [
    #'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root'
    'file:/scratchdisk2/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root'
    ]

process = makeZprime2muAnalysisProcess(files)
attachAnalysis(process, 'Zprime2muExample')

process.Zprime2muExample.nBins     = cms.int32(40)
process.Zprime2muExample.lowerMass = cms.double(400)
process.Zprime2muExample.upperMass = cms.double(2000)
