import FWCore.ParameterSet.Config as cms

# Hack to set the python path for now in 1_6_X.
import sys, os
sys.path.insert(0, os.path.join(os.environ['CMSSW_BASE'],
                                'src/SUSYBSMAnalysis/Zprime2muAnalysis/data'))

from Zprime2muAnalysisCommon_cff import *
from Zprime2muResolution_cfi import *

files = [
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_100_11110001.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_100_22220001.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_100_33330001.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_11110020.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_11110300.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_11114000.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_22220020.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_22220300.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_22224000.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_33330020.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_33330300.root',
    'rfio:/castor/cern.ch/user/p/pivarski/Zprime_206_FEVTSIM/Zprime1TeV_400_33334000.root'
    ]

process = makeZprime2muAnalysisProcess(files,
                                       useHEEPSelector=False,
                                       useTrigger=False,
                                       useOtherMuonRecos=False,
                                       recoverBrem=False,
                                       disableElectrons=True)
#   , doingElectrons=True,
#   flavorsForDileptons=diElectrons)

attachResolution(process)

process.Zprime2muResolution.verbosity = 2
process.Zprime2muResolution.dateHistograms = False

# Example of how to replace one of the parameter sets:
'''
process.Zprime2muResolution.exampleSet = cms.PSet(
    peakMass = cms.double(200),
    lowerMassWin = cms.double(0),
    upperMassWin = cms.double(300),
    binSize = cms.int32(25),
    maxTrigMass = cms.double(0.3))
process.Zprime2muResolution.dataSet = 'exampleSet'
'''
