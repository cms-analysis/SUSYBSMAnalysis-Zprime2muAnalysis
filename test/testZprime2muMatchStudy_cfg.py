import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root'
    ]

process = makeZprime2muAnalysisProcess(files)
attachAnalysis(process, 'Zprime2muMatchStudy', isRecLevelAnalysis=True)
