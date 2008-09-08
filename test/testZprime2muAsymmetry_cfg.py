import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAsymmetry_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/dy_above400_RECO.root'
    ]

process = makeZprime2muAnalysisProcess(files)
attachAsymmetry(process)

#process.Zprime2muAsymmetry.numFits = 4
#process.Zprime2muAsymmetry.useCachedParams = true
