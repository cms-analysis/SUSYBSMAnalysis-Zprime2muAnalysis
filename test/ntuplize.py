import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *

files = ['file:/scratchdisk3/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_1.root']
process = makeZprime2muAnalysisProcess(files)
dumpNtuple(process, files[0].replace('.root', '.zp2muntuple.root'))
