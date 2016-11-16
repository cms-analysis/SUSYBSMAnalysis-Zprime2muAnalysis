#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring('file:PlaceHolder.root'))

# Load services needed to run the PAT.
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = cms.string('PlaceHolder::All')

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(limit = cms.untracked.int32(-1)) 

## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)

process.load('PhysicsTools.PatAlgos.patSequences_cff')

# Define the output file with the output commands defining the
# branches we want to have in our PAT tuple.
process.out = cms.OutputModule('PoolOutputModule',
                               fileName = cms.untracked.string('pat.root'),
                               # fileName = cms.untracked.string('file:PlaceHolder.root'),
                               SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep patElectrons_*_*_*',
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
                                                                      'keep *_cleanPatMuonsTriggerMatch_*_*', 
                                                                      'keep *_patTrigger_*_*', 
                                                                      'keep *_patTriggerEvent_*_*', 
                                                                      'keep *_patMETsPF_*_*',
                                                                      )
                               )

# PAT taus
del process.patTaus.tauIDSources.againstElectronMVA6Raw
#del process.patTaus.tauIDSources.againstMuonMedium
process.cleanPatTaus.preselection = process.cleanPatTaus.preselection.value().replace('againstElectronVLooseMVA5', 'againstElectronLooseMVA6')
process.cleanPatTaus.preselection = process.cleanPatTaus.preselection.value().replace('againstMuonMedium', 'againstMuonTight')

# PAT muons 
process.patMuons.embedTrack = True
process.selectedPatMuons.cut = "isTrackerMuon || isGlobalMuon"
process.countPatMuons.minNumber = cms.uint32(1)

# PAT trigger info
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi')
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger, switchOnTriggerMatchEmbedding
switchOnTrigger( process )
switchOnTriggerMatchEmbedding( process,
                                triggerProducer = 'patTrigger', # this is already the default setting
                                triggerMatchers = [ 'muonTriggerMatchHLTMuons' ]
                               )
# PAT Met and Jets
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection 
addMETCollection(process, labelName='patMETsPF', metSource='pfMetT1')

from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection #to be checked
switchJetCollection(process, 
                    jetSource = cms.InputTag('ak4PFJets'),
                    jetCorrections = ('AK4PF', cms.vstring(['L1FastJet', 'L2Relative','L3Absolute']), 'Type-1'),
                    btagDiscriminators = ['jetBProbabilityBJetTags',
                                          'jetProbabilityBJetTags',
                                          'trackCountingHighPurBJetTags',
                                          'trackCountingHighEffBJetTags',
                                          'simpleSecondaryVertexHighEffBJetTags',
                                          'simpleSecondaryVertexHighPurBJetTags',
                                          'combinedSecondaryVertexV2BJetTags'],
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

# Met filters
process.load("PhysicsTools.PatAlgos.slimming.metFilterPaths_cff")
process.goodDataHBHENoiseFilter       = cms.Path(process.HBHENoiseFilter)
process.goodDataHBHEIsoNoiseFilter       = cms.Path(process.HBHENoiseIsoFilter)
process.goodDataCSCTightHaloFilter    = cms.Path(process.globalTightHalo2016Filter)
process.goodDataEcalTPFilter          = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.goodDataEeBadScFilter         = cms.Path(process.eeBadScFilter)

##if you want to filter the event: define a sequence and include it in the pocess.p
#process.goodDataMETFilter =  cms.Sequence(process.HBHENoiseFilter * process.CSCTightHaloFilter * process.hcalLaserEventFilter * process.EcalDeadCellTriggerPrimitiveFilter * process.trackingFailureFilter * process.eeBadScFilter * process.ecalLaserCorrFilter *process.trkPOGFilters)

##if you want just to tag the event: define a path
process.goodDataMETFilter =  cms.Path(process.HBHENoiseFilter *
				      process.HBHENoiseIsoFilter *  
                                      process.globalTightHalo2016Filter * 
                                      process.EcalDeadCellTriggerPrimitiveFilter * 
                                      process.eeBadScFilter)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
process.goodDataHLTPhysicsDeclared = cms.Path(process.hltPhysicsDeclared)
process.goodDataPrimaryVertexFilter = cms.Path(process.primaryVertex)
process.goodDataNoScraping = cms.Path(process.noscraping)
process.goodDataAll = cms.Path(process.hltPhysicsDeclared * process.primaryVertex) # * process.noscraping)

process.outpath = cms.EndPath(process.out) 
#print process.dumpPython()
