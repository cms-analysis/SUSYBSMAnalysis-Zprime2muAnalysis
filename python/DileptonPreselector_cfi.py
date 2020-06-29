import FWCore.ParameterSet.Config as cms

dileptonPreselector = cms.EDFilter("DileptonPreselector",
	muons = cms.InputTag("slimmedMuons"),
	nMuons = cms.double(2),
    # This was 40 GeV for some unknown reason. It needs to be low 
    # enough to be able to do the Z0 normalization measurement. 
    # That is, some value less than whatever prescaled trigger
    # we use (HLT_Mu27). 
	ptCut = cms.double(20),	 

)
