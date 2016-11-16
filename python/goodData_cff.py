import FWCore.ParameterSet.Config as cms

# The below is mostly copied from
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes ; the
# user is responsible for checking that what's used is up-to-date and
# appropriate to their analysis.

hltPhysicsDeclared = cms.EDFilter('HLTPhysicsDeclared',
                                  invert = cms.bool(False),
                                  L1GtReadoutRecordTag = cms.InputTag('gtDigis')
                                  )

primaryVertex = cms.EDFilter('GoodVertexFilter',
                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                   minimumNDOF = cms.uint32(4),
                                   maxAbsZ = cms.double(24),
                                   maxd0 = cms.double(2)
                                   )

noscraping = cms.EDFilter('FilterOutScraping',
                          applyfilter = cms.untracked.bool(True),
                          debugOn = cms.untracked.bool(False),
                          numtrack = cms.untracked.uint32(10),
                          thresh = cms.untracked.double(0.25)
                          )
from METFilterMiniAOD_cfi import defaultSelector

primaryVertexMiniAOD = defaultSelector.clone()
primaryVertexMiniAOD.flag = "Flag_goodVertices"

beamHaloMiniAOD = defaultSelector.clone()
beamHaloMiniAOD.flag = "Flag_globalTightHalo2016Filter"

hbheNoiseMiniAOD = defaultSelector.clone()
hbheNoiseMiniAOD.flag = "Flag_HBHENoiseFilter"

hbheIsoNoiseMiniAOD = defaultSelector.clone()
hbheIsoNoiseMiniAOD.flag = "Flag_HBHENoiseIsoFilter"

ecalDeadCellPrimitiveMiniAOD = defaultSelector.clone()
ecalDeadCellPrimitiveMiniAOD.flag = "Flag_EcalDeadCellTriggerPrimitiveFilter"

eeBadScMiniAOD = defaultSelector.clone()
eeBadScMiniAOD.flag = "Flag_eeBadScFilter"

metFilters = [beamHaloMiniAOD, hbheNoiseMiniAOD, hbheIsoNoiseMiniAOD, ecalDeadCellPrimitiveMiniAOD, eeBadScMiniAOD]

def addNewFilters(process,filterList):

	process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
	process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
	process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

	process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
	process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
	process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

	filterList.append(process.BadPFMuonFilter)
	filterList.append(process.BadChargedCandidateFilter)

	return filterList
