import FWCore.ParameterSet.Config as cms

# Hack to set the python path for now in 1_6_X.
import sys, os
sys.path.insert(0, os.path.join(os.environ['CMSSW_BASE'],
                                'src/SUSYBSMAnalysis/Zprime2muAnalysis/data'))

from Zprime2muAnalysisCommon_cff import *
from Zprime2muAnalysis_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root'
    ]

process = makeZprime2muAnalysisProcess(files)
#   ,doingElectrons=True,
#   flavorsForDileptons=diElectrons)

attachAnalysis(process)

#print process.dumpConfig()
