import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muHistos_cfi import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muResolution_cfi import *

cast_dir = 'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2'
mass, lcut = 1000, 400
#mass, lcut = 1300, 600
sample = 'ZPSSMmumu_M%(mass)i_Mcut%(lcut)i' % locals()
gtag = "MC_31X_V3"
#gtag = "STARTUP31X_V2_EX"
num_files = 25

files = ['%(cast_dir)s/%(sample)s-%(gtag)s/PYTHIA6_%(sample)s_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_%(gtag)s_%(i)i.root' % locals()
         for i in xrange(1, num_files + 1)]

process = makeZprime2muAnalysisProcess(files, maxEvents=1000) #, useTrigger=False)

attachHistos(process) # Run Histos, too.
attachResolution(process, verbosity=2)

#for module in [process.Zprime2muHistos, process.Zprime2muResolution]:
    #setattr(module, 'leptonsFromDileptons', True)
    #setattr(module, 'dataSet', 'DY40')

#print process.dumpConfig()
