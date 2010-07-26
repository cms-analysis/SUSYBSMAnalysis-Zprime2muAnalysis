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


## ## Define some utilities to declare bins easily
## def nBins(n,min,max): return cms.vdouble(*[min + (max-min)/n*i for i in range(n+1)])
## def evenBins(min,max,delta): 
##     ret = cms.vdouble(min)
##     x = min
##     while x < max - 1e-4: # need a small hint otherwise for some numbers it will overstep due to numerical resolution
##         x += delta
##         ret.append(x)
##     return ret 

## rebinFactor=1.0
## process.ggPlots = cms.EDAnalyzer("Zprime2muAnalysisPlots",
##     #muons     = cms.InputTag('gMuons'),
##     dilepton_src = cms.InputTag('zToMuMuGG'),
##     selection = cms.string(""), #"isTrackerMuon && muonID('TMLastStationAngTight')"),
##     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
##     #normalization   = cms.InputTag("none"), # optional!
##     # ---- Kinematics ----
##     massBins = nBins( 100, 0., 100.),
##     ptBins = evenBins( 0, 15, 0.25),
##     pBins  = evenBins( 0, 30, 0.5),
##     etaBins = evenBins( -2.6, 2.6, 0.2),
##     phiBins = evenBins(-3.2,  3.2, 0.2),
##     chargeBins = cms.vdouble(-2,0,2),
##     # ---- Vertex ----
##     dxyFineBins = evenBins(-0.5,0.5, 0.005), #  50um
##     dzFineBins  = evenBins(-1.0,1.0, 0.010), # 100um
##     dxyCoarseBins = evenBins(-10,10,0.1), # 1mm
##     dzCoarseBins  = evenBins(-30,30,0.1), # 1mm
##     # ---- Tracks ----
##     pixelHitsBins   = nBins(8,0,8),
##     pixelLayersBins     = nBins(5,0,5),                         
##     trackerHitsBins = nBins(33,0,33),
##         trackerLostHitsBins = nBins(10,0,10),
##     muonHitsBins    = nBins(25,0,50),
##     muonStationHitsBins = nBins(20,0,20),
##     muonBadHitsBins     = nBins(20,0,20),
##     globalHitsBins  = nBins(40,0,80),
##     trackerChi2nBins = evenBins(0, 10, 0.2),
##     muonChi2nBins    = evenBins(0, 10, 0.2),
##     globalChi2nBins  = evenBins(0, 10, 0.2),
##     # ---- Isolation ----
##     isolationBins = evenBins(0, 5, .5),
##                                    relIsoBins    = evenBins(0, .5, .025 * rebinFactor),
##         # ---- Muon ID ----
##         muonStationsBins    = nBins(5,0,5), 
##         segmentMatchesBins = nBins(12,0,12),
##         segmentCompatBins  = evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
##         caloCompatBins     = evenBins(0, 1 + 0.1*rebinFactor, 0.1 * rebinFactor), # need one bin for ">= 1.0"
##         boolBins = nBins(2,-0.5,1.5),
##         zBins = nBins(100,-500.,500.),
##         rBins = nBins(100,0.,500.),
##         rzXBins = cms.uint32(1000),
##         rzXRange = cms.vdouble(-500.,500.),
##         rzYBins = cms.uint32(500),
##         rzYRange = cms.vdouble(0.,500.),

## )

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
     process.ttPlots)
    )
