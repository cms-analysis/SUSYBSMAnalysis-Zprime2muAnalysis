import FWCore.ParameterSet.Config as cms

EventCounter = cms.EDAnalyzer('EventCounter',
	genInfoTag = cms.InputTag("generator"),
)
