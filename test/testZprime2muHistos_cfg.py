import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muHistos_cfi import *

files = [
    'file:/scratchdisk3/tucker/ZPSSMmumu_M1000_Mcut400_10TeV_IDEAL_V9_RAW2DIGI_RECO_1-10.root'
    ]

process = makeZprime2muAnalysisProcess(files, skipPAT=True, disableElectrons=True)
attachHistos(process, verbosity=2) #, leptonsFromDileptons=True)
