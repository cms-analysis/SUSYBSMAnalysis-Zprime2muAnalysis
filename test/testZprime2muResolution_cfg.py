import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muHistos_cfi import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muResolution_cfi import *

mass = 1000
lcut = 400
gtag = "MC_31X_V3"
#gtag = "STARTUP31X_V2_EX"

dir = "rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/"
sample = "ZPSSMmumu_M" + str(mass) + "_Mcut" + str(lcut) + "-" + gtag + "/"
attr = "_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_"
file = "PYTHIA6_ZPSSMmumu_M" + str(mass) + "_Mcut" + str(lcut) + attr + gtag

files = [
    dir + sample + file + "_1.root",
    dir + sample + file + "_2.root",
    dir + sample + file + "_3.root",
    dir + sample + file + "_4.root",
    dir + sample + file + "_5.root",
    dir + sample + file + "_6.root",
    dir + sample + file + "_7.root",
    dir + sample + file + "_8.root",
    dir + sample + file + "_9.root",
    dir + sample + file + "_10.root",
    dir + sample + file + "_11.root",
    dir + sample + file + "_12.root",
    dir + sample + file + "_13.root",
    dir + sample + file + "_14.root",
    dir + sample + file + "_15.root",
    dir + sample + file + "_16.root",
    dir + sample + file + "_17.root",
    dir + sample + file + "_18.root",
    dir + sample + file + "_19.root",
    dir + sample + file + "_20.root",
    dir + sample + file + "_21.root",
    dir + sample + file + "_22.root",
    dir + sample + file + "_23.root",
    dir + sample + file + "_24.root",
    dir + sample + file + "_25.root"
    ]

process = makeZprime2muAnalysisProcess(files, skipPAT=True, disableElectrons=True, maxEvents=1000) #, useTrigger=False)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

attachHistos(process) # Run Histos, too.
attachResolution(process, verbosity=2)

#for module in [process.Zprime2muHistos, process.Zprime2muResolution]:
    #setattr(module, 'leptonsFromDileptons', True)
    #setattr(module, 'dataSet', 'DY40')

#print process.dumpConfig()
