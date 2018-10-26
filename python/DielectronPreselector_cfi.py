import FWCore.ParameterSet.Config as cms

dielectronPreseletor = cms.EDFilter("DielectronPreselector",
	electrons = cms.InputTag("slimmedElectrons"),
	nElectrons = cms.double(2),
	ptCut = cms.double(20),	

)
