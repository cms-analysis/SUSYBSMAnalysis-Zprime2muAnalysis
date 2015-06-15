#!/usr/bin/env python

# Parameters you usually want to change are here at the top of the
# file.

# Change the following to False to run on data.
input_is_MC = True

# Change the following parameters appropriately for the samples over
# which you wish to run, most especially the globalTag (since the list
# of files will be configured by CRAB itself so it doesn't matter what
# gets put here for it).
globalTag = 'START53_V11::All'
files = ['file:/eos/uscms/store/mc/Summer12_DR53X/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/AODSIM/PU_S10_START53_V7C1-v1/10000/A010B486-703D-E211-8903-78E7D1E4B6E8.root']
#files = ['root://cmsxrootd-site.fnal.gov/store/mc/Summer12_DR53X/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/AODSIM/PU_S10_START53_V7C1-v1/10000/A010B486-703D-E211-8903-78E7D1E4B6E8.root']
maxEvents = -1

########################################################################

# Load and configure the process appropriately.

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.GlobalTag.globaltag = globalTag
process.source.fileNames = files
process.maxEvents.input = maxEvents

# Do the right thing depending on whether MC truth is supposed to be
# available.
if not input_is_MC:
    from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
    removeMCUse(process)

# Our path. The name 'p' is expected by the configuration of the
# OutputModule in PATTuple_cfg.
process.p = cms.Path(process.type0PFMEtCorrection * process.patDefaultSequence)

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import AODOnly
AODOnly(process)

########################################################################

# Changing things/uncommenting below here can be done by slightly more
# advanced users :-) You can use a bunch of core tools of PAT and some
# Exotica tools to tailor your PAT configuration; for a few, examples
# are commented out in the next several lines.

#from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD, removeAllPATObjectsBut, removeSpecificPATObjects
#restrictInputToAOD(process)
#removeAllPATObjectsBut(process, ['Muons', 'Electrons'])
#removeSpecificPATObjects(process, ['Taus'])

# Some MC samples have the HLT process name different from "HLT".
#from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName
#switchHLTProcessName(process, "REDIGI311X")

# You can change the default muon cut from "isGlobalMuon ||
# isTrackerMuon" to your taste. Some examples (only uncomment one!):
#process.selectedPatMuons.cut = 'isTrackerMuon && pt > 1 && p > 2.5 && innerTrack.hitPattern.numberOfValidTrackerHits > 12 && innerTrack.normalizedChi2 < 5 && abs(dB) < 0.5 && abs(dZ) < 5 && muonID("TMLastStationAngTight")'
#process.selectedPatMuons.cut = 'muonID("GlobalMuonPromptTight") && muonID("TMOneStationLoose") && (globalTrack.hitPattern.numberOfValidMuonCSCHits + globalTrack.hitPattern.numberOfValidMuonDTHits) >= 1 && innerTrack.hitPattern.trackerLayersWithMeasurement >= 6'

# You may only want to keep electrons that pass the HEEP cuts, but
# does this screw up the PAT jet cleaning? Needs to be studied.
#process.selectedPatElectrons.cut = 'userInt("HEEPId") == 0'

# If you're studying e.g. e/mu dileptons, you may want to change the
# requirement from one/two muons to two leptons (= electrons + muons).
#process.countPatMuons.minNumber = 0
#process.countPatLeptons.minNumber = 2
