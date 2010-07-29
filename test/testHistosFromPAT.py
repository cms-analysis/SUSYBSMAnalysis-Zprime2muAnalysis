#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
process.leptons.muon_cuts = ''

process.gMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("leptons:muons"),
    cut = cms.string("isGlobalMuon"), 
)

process.tMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("leptons:muons"),
    cut = cms.string("isTrackerMuon && !isGlobalMuon"), 
)

process.zToMuMuGG = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('gMuons@+ gMuons@-'),
    cut = cms.string('0.0 < mass < 20000.0'),
    name = cms.string('zToMuMuGG'),
    roles = cms.vstring('muon1', 'muon2')
)

process.zToMuMuGT = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('gMuons@+ tMuons@-'),
    cut = cms.string('0.0 < mass < 20000.0'),
    name = cms.string('zToMuMuGT'),
    roles = cms.vstring('muon1', 'muon2')
)

process.zToMuMuTT = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tMuons@+ tMuons@-'),
    cut = cms.string('0.0 < mass < 20000.0'),
    name = cms.string('zToMuMuTT'),
    roles = cms.vstring('muon1', 'muon2')
)



from SUSYBSMAnalysis.Zprime2muAnalysis.inclusiveDiMuonPlots_cfi import makeInclusiveDiMuonPlots
commonInputs = cms.PSet(
    dilepton_src = cms.InputTag('zToMuMuGG'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    selection = cms.string(""),
)
process.ggPlots = cms.EDAnalyzer("Zprime2muAnalysisPlots",
                                 makeInclusiveDiMuonPlots(),
                                 commonInputs)
process.gtPlots = process.ggPlots.clone(dilepton_src = "zToMuMuGT")
process.ttPlots = process.ggPlots.clone(dilepton_src = "zToMuMuTT")


from UserCode.Examples.inclusiveMuonPlotsMRTU_cfi import makeInclusiveMuonPlots;
commonInputsLepton = cms.PSet(
    muons     = cms.InputTag('cleanPatMuons'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
)
process.trackerMuons = cms.EDAnalyzer("InclusiveMuonPlotsMRTU",
    makeInclusiveMuonPlots(),
    commonInputsLepton,
    selection = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
)
process.globalMuons = process.trackerMuons.clone(
    selection = "isGlobalMuon"
    )
process.standAloneMuons = process.trackerMuons.clone(
    selection = "isStandAloneMuon"
    )


process.p = cms.Path(
#    process.hltFilter +
    process.Zprime2muAnalysisSequence
    (process.gMuons *
     process.tMuons ) *
    (process.zToMuMuGG *
    process.zToMuMuGT *
    process.zToMuMuTT) *
    (process.ggPlots +
     process.gtPlots +
     process.ttPlots) *
    (process.trackerMuons *
     process.globalMuons *
     process.standAloneMuons 
     )
    )
