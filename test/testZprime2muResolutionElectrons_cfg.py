import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muResolution_cfi import *

files = [
    '/store/mc/2007/12/11/CSA07-ZPSSMee_M1000_Mcut400-1197356198/0025/06E3882D-DAA7-DC11-B45E-001617C3B716.root'
    ]

process = makeZprime2muAnalysisProcess(files, doingElectrons=True, flavorsForDileptons=diElectrons)
attachResolution(process)

process.Zprime2muResolution.verbosity = 2

