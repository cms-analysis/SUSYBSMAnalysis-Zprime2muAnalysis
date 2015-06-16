#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:PlaceHolder.root'))

# Load services needed to run the PAT.
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag 
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup') 
process.GlobalTag.globaltag = cms.string('PlaceHolder::All')

# Configure the MessageLogger ~sanely. Also direct it to let the PAT
# summary tables be reported -- nice to see how many events had no
# muons, how many had no "selected"/"clean" muons, etc.
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(limit = cms.untracked.int32(-1))

## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)

# Load the PAT modules and sequences, and configure them as we
# need. See the individual functions for their documentation.  MC use
# is assumed by default, and should be removed after everything's
# configured in the top-level config using removeMCUse() tool if
# running on data.  (This is due to the design of the PAT: easier to
# do it in this order rather than adding things for MC use later.)
process.load('PhysicsTools.PatAlgos.patSequences_cff')

# Define the output file with the output commands defining the
# branches we want to have in our PAT tuple.
process.out = cms.OutputModule('PoolOutputModule',
                               fileName = cms.untracked.string('pat.root'),
                               # If your path in your top-level config is not called 'p', you'll need
                               # to change the next line. In test/PATTuple.py, 'p' is used.
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('p')),

                               outputCommands = cms.untracked.vstring(
                                   'drop *',
                                   'keep patElectrons_cleanPatElectrons__*',
                                   'keep patMuons_cleanPatMuons__*',
                                   'keep patJets_cleanPatJets__*',
                                   'keep patPhotons_cleanPatPhotons__*',
                                   'keep patMETs_patMETs*__PAT',
                                   'keep recoGenParticles_prunedMCLeptons_*_*',
                                   'keep recoGenJets_selectedPatJets_genJets_*',
                                   'keep recoBaseTagInfosOwned_selectedPatJets_tagInfos_*',
                                   'keep GenEventInfoProduct_*_*_*',
                                   'keep GenRunInfoProduct_*_*_*',
                                   'keep *_offlineBeamSpot_*_*',
                                   'keep *_offlinePrimaryVertices_*_*',
                                   'keep edmTriggerResults_TriggerResults__HLT*',
                                   'keep edmTriggerResults_TriggerResults__REDIGI*',
                                   'keep L1GlobalTriggerObjectMaps_l1L1GtObjectMap_*_*', # drop later if embedding of L1 into PAT works fine
                                   'keep L1GlobalTriggerReadoutRecord_gtDigis__RECO',
                                   'keep *_hltTriggerSummaryAOD__HLT*',
                                   'keep *_hltTriggerSummaryAOD__REDIGI*',
                                   'keep edmTriggerResults_TriggerResults__PAT', # for post-tuple filtering on the goodData paths
                                   'keep PileupSummaryInfos_addPileupInfo_*_*',   # may be needed for pile-up reweighting
                                                                      #
                                   'keep *_cleanPatMuonsTriggerMatch_*_*',
                                   'keep *_patTrigger_*_*', # keep these two for now, for Level-1 decisions
                                   'keep *_patTriggerEvent_*_*',
                                   )
                               )

# PAT electrons
from PATTools import addHEEPId
addHEEPId(process)

# PAT taus
del process.patTaus.tauIDSources.againstElectronMVA5raw
del process.patTaus.tauIDSources.againstMuonMedium
process.cleanPatTaus.preselection = process.cleanPatTaus.preselection.value().replace('againstMuonMedium', 'againstMuonTight') #now Tight is the default choice

# PAT muons
# Embed the tracker tracks (by default, every other track is already
# embedded).
process.patMuons.embedTrack = True
# Follow VBTF for now in using the beamspot for "correcting" dxy,
# instead of the primary vertex.
process.patMuons.usePV = False
# Define simple quality cuts for muons (analysis cuts to be done at the analysis level).

process.selectedPatMuons.cut = 'isGlobalMuon || isTrackerMuon'

# Filter out events with no selected muons. (countPatMuons counts
# those muons in the cleanPatMuons collection, which by default is
# just a copy of the selectedPatMuons collection.) If, e.g., wanting
# to use the PAT tuple for efficiency studies, or running on just
# electrons, need to change appropriately in your top-level config
# (perhaps using countPatLeptons instead).
process.countPatMuons.minNumber = 1

# Add MET and jets.
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
addMETCollection(process, labelName='patMETsPF', metSource='pfMetT1')

from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
switchJetCollection(process,
                    jetSource = cms.InputTag('ak5PFJets'),
                    jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative','L3Absolute']), 'Type-1'),
                    btagDiscriminators = ['jetBProbabilityBJetTags',
                                          'jetProbabilityBJetTags',
                                          'trackCountingHighPurBJetTags',
                                          'trackCountingHighEffBJetTags',
                                          'simpleSecondaryVertexHighEffBJetTags',
                                          'simpleSecondaryVertexHighPurBJetTags',
                                          'combinedSecondaryVertexBJetTags',
                                          ],
                    getJetMCFlavour = False,
                    )

# Make a collection of muons with our final selection applied so that
# the muon-jet cleaning will use only muons passing those cuts. This
# muon collection is not saved in the output.
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff import loose_cut
process.muonsForJetCleaning = process.selectedPatMuons.clone(cut = loose_cut.replace('pt > 45', 'pt > 30'))
process.patDefaultSequence.replace(process.selectedPatMuons, process.selectedPatMuons * process.muonsForJetCleaning)
process.cleanPatJets.checkOverlaps.muons.src = 'muonsForJetCleaning'
process.cleanPatJets.checkOverlaps.muons.deltaR = 0.2
process.cleanPatJets.checkOverlaps.muons.requireNoOverlaps = True
process.cleanPatJets.finalCut = 'pt > 30.0'

# PAT Trigger info
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger, switchOnTriggerMatchEmbedding
switchOnTrigger( process )
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi')
switchOnTriggerMatchEmbedding(process,
                              #triggerProducer = 'patTrigger',
                              triggerMatchers=['muonTriggerMatchHLTMuons']
                              )


# Instead of filtering out events at PAT-tupling time based on things
# like GoodVertex and NoScraping, schedule separate paths for all the
# "good data" filters so that the results of them get stored in a
# small TriggerResults::PAT object. This can be read and used to
# filter events in the analyzer process using e.g. the filter in
# Zprime2muAnalysis_cff.py. (This is useful so we don't have to keep
# around the entire generalTracks collection to run the NoScraping
# filter later, for example.)
#
# Make one path for each (a very small storage burden) so they can be
# accessed separately in the TriggerResults object;

# MET filters.
#process.load("RecoMET.METFilters.metFilters_cff") #this gives an error
process.load("PhysicsTools.PatAlgos.slimming.metFilterPaths_cff")
process.goodDataTrackingFailureFilter = cms.Path(process.trackingFailureFilter)
process.trackingFailureFilter.VertexSource = cms.InputTag('offlinePrimaryVertices')
process.goodDataTrackingPOGFilter     = cms.Path(process.trkPOGFilters)
process.goodDataCSCTightHaloFilter    = cms.Path(process.CSCTightHaloFilter)
process.goodDataEcalTPFilter          = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.goodDataEeBadScFilter         = cms.Path(process.eeBadScFilter)
process.goodDataEcalLaserFilter       = cms.Path(process.ecalLaserCorrFilter)
process.goodDataHBHENoiseFilter       = cms.Path(process.HBHENoiseFilter)
process.goodDataHcalLaserFilter       = cms.Path(process.hcalLaserEventFilter)

##if you want to filter the event: define a sequence and include it in the pocess.p
#process.goodDataMETFilter =  cms.Sequence(process.HBHENoiseFilter * process.CSCTightHaloFilter * process.hcalLaserEventFilter * process.EcalDeadCellTriggerPrimitiveFilter * process.trackingFailureFilter * process.eeBadScFilter * process.ecalLaserCorrFilter *process.trkPOGFilters)
##if you want just to tag the event: define a path
process.goodDataMETFilter =  cms.Path(process.HBHENoiseFilter * process.CSCTightHaloFilter * process.hcalLaserEventFilter * process.EcalDeadCellTriggerPrimitiveFilter * process.trackingFailureFilter * process.eeBadScFilter * process.ecalLaserCorrFilter *process.trkPOGFilters)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
process.goodDataHLTPhysicsDeclared = cms.Path(process.hltPhysicsDeclared)
process.goodDataPrimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.goodDataNoScraping = cms.Path(process.noscraping)
# The "All" path isn't necessary because it could be emulated using
# the AND of all of the separate ones, but it's nice for convenience.
# Only include primary vertex and HLT declared for now.
process.goodDataAll = cms.Path(process.hltPhysicsDeclared * process.primaryVertexFilter)#* process.noscraping)

process.outpath = cms.EndPath(process.out)
#print process.dumpPython()
