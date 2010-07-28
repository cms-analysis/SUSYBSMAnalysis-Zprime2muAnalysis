#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

process = cms.Process('Zprime2muAnalysis')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:patTuple.root'))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service('TFileService', fileName=cms.string('zp2mu_histos.root'))

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.HLTPaths = ['HLT_Mu9', 'HLT_DoubleMu3']
process.hltFilter.andOr = True # == OR

process.leptons = cms.EDProducer('Zprime2muLeptonProducer',
                         muon_src = cms.InputTag('cleanPatMuons'),
                         electron_src = cms.InputTag('cleanPatElectrons'),
                         muon_cuts = cms.string(''), #'pt > 20. && isolationR03.sumPt < 10'),
                         electron_cuts = cms.string(''),
                         muon_track_for_momentum = cms.string('pmc'),
                         muon_photon_match_src = cms.InputTag('muonPhotonMatch')
                         )

process.load("SUSYBSMAnalysis.Zprime2muAnalysis.MuonPhotonMatch_cff")
#from SUSYBSMAnalysis.Zprime2muAnalysis.adam_photonMatch_cff import addUserData as addPhotonMatch
#addPhotonMatch(process.leptons)


process.gMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("leptons:muons"),
    cut = cms.string("isGlobalMuon "), 
)

process.tMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("leptons:muons"),
    cut = cms.string("isTrackerMuon && !isGlobalMuon"), 
)

process.zToMuMuGG = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('gMuons@+ gMuons@-'),
    #decay = cms.string('cleanPatMuons@+ cleanPatMuons@-'),
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



from SUSYBSMAnalysis.Zprime2muAnalysis.inclusiveDiMuonPlots_cfi import makeInclusiveDiMuonPlots;
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
    process.muonPhotonMatch *
    process.leptons *
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
