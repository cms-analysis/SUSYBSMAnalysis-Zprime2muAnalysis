#!/usr/bin/env python

# Parameters you usually want to change are here at the top of the
# file.

# Change the following to False to run on data.
input_is_MC = True

# Change the following parameters appropriately for the samples over
# which you wish to run, most especially the globalTag (since the list
# of files will be configured by CRAB itself so it doesn't matter what
# gets put here for it).
globalTag = 'START36_V10::All'
files = ['/store/mc/Summer10/Zmumu/GEN-SIM-RECO/START36_V9_S09-v1/0045/FC750941-147B-DF11-BA8C-00261834B569.root']

########################################################################

# Load and configure the process appropriately.

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.GlobalTag.globaltag = globalTag
process.source.fileNames = files

# You can make some typical event selection cuts such as good vertex,
# beam scraping, and L1 technical bits by including it in the path.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')

# Do the right thing depending on whether MC truth is supposed to be
# available.
if input_is_MC:
    # Modify L1 trigger bits if using MC (bit 0 not defined for MC).
    process.L1T1.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
else:
    from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
    removeMCUse(process)

# Our path. The name 'p' is expected by the configuration of the
# OutputModule in PATTuple_cfg.
process.p = cms.Path(process.goodData * process.patDefaultSequence)

########################################################################

# Changing things/uncommenting below here can be done by slightly more
# advanced users :-)

# You can use a bunch of core tools of PAT and some Exotica tools to
# tailor your PAT configuration; for a few examples are commented out
# in the next several lines.
#from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD, removeAllPATObjectsBut, removeSpecificPATObjects
#restrictInputToAOD(process)
#removeAllPATObjectsBut(process, ['Muons', 'Electrons'])
#removeSpecificPATObjects(process, ['Taus'])

# You can still run PAT in the 36X version on input files produced
# within the 35X series. This requires some reconfigurations that are
# done by the PAT for us with a canned function. Uncomment the next
# two lines to use it.
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run36xOn35xInput
#run36xOn35xInput(process)

# You can change the default muon cut from "isGlobalMuon ||
# isTrackerMuon" to your taste.
#process.selectedPatMuons.cut = 'isTrackerMuon && pt > 1 && p > 2.5 && innerTrack.hitPattern.numberOfValidTrackerHits > 12 && innerTrack.normalizedChi2 < 5 && abs(dB) < 0.5 && abs(dZ) < 5 && muonID("TMLastStationAngTight")'
#process.selectedPatMuons.cut = 'muonID("GlobalMuonPromptTight") && muonID("TMOneStationLoose") && (globalTrack.hitPattern.numberOfValidMuonCSCHits + globalTrack.hitPattern.numberOfValidMuonDTHits) >= 1 && innerTrack.hitPattern.trackerLayersWithMeasurement >= 6'

# You may only want to bother with electrons that pass the HEEP cuts,
# but will this screw up the PAT jet cleaning? Needs to be studied.
#process.selectedPatElectrons.cut = 'userInt("HEEPId") == 0'

# If you're studying e.g. e/mu dileptons, you may want to 
#process.countPatMuons.minNumber = 0
#process.countPatLeptons.minNumber = 2

# Options for controlling how CMSSW works.
#process.source.noEventSort = cms.untracked.bool(True)
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
#process.maxEvents.input = 1000
#process.out.fileName = 'myTuple.root'
#process.options.wantSummary = False # to suppress the long output at the end of the job

