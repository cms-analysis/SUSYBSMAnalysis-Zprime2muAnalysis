import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muHistos_cfi import *

files = [
    'file:/scratchdisk3/tucker/ZPSSMmumu_M1000_Mcut400_10TeV_IDEAL_V9_RAW2DIGI_RECO_1-10.root'
#    'file:/scratchdisk3/tucker/RelValZEE-2_1_10.root'
    ]

process = makeZprime2muAnalysisProcess(files, 5, skipPAT=True, electrons=([None]*9))
attachData(process)

process.Zprime2muHistos.verbosity = 2
process.Zprime2muHistos.dateHistograms = False
